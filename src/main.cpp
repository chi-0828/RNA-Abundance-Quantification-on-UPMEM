#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>
#include <algorithm>
#include <limits>

#include <cstdio>

#include <zlib.h>

#include "common.h"
#ifdef DPU_HEAD
#define DPU_HEAD 0
#else
#define DPU_HEAD 1
#endif
#include "ProcessReads.h"
#include "KmerIndex.h"
#include "Kmer.hpp"
#include "MinCollector.h"
#include "EMAlgorithm.h"
#include "weights.h"
#include "Inspect.h"
#include "Bootstrap.h"
#include "H5Writer.h"
#include "PlaintextWriter.h"
#include "GeneModel.h"
#include "Merge.h"


//#define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"

using namespace std;


int my_mkdir(const char *path, mode_t mode) {
  #ifdef _WIN64
  return mkdir(path);
  #else
  return mkdir(path,mode);
  #endif
}

bool checkFileExists(std::string fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}


void ParseOptionsPseudo(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int single_flag = 0;
  int strand_flag = 0;
  int strand_FR_flag = 0;
  int strand_RF_flag = 0;
  int pbam_flag = 0;
  int gbam_flag = 0;
  int umi_flag = 0;
  int quant_flag = 0;
  int bus_flag = 0;
  int dpu_flag = 64;
  int dpurk_flag = 1;

  const char *opt_string = "t:i:l:s:o:b:d:r:u:g:n";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"single", no_argument, &single_flag, 1},
    //{"strand-specific", no_argument, &strand_flag, 1},
    {"fr-stranded", no_argument, &strand_FR_flag, 1},
    {"rf-stranded", no_argument, &strand_RF_flag, 1},
    {"pseudobam", no_argument, &pbam_flag, 1},
    {"quant", no_argument, &quant_flag, 1},
    {"bus", no_argument, &bus_flag, 1},
    {"num", no_argument, 0, 'n'},
    {"umi", no_argument, &umi_flag, 'u'},
    {"batch", required_argument, 0, 'b'},
    {"dpu", required_argument, &dpu_flag, 'd'},
    {"rank", required_argument, &dpurk_flag, 'r'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"gtf", required_argument, 0, 'g'},
    {"index", required_argument, 0, 'i'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'd': {
      opt.dpu = atoi(optarg);
      break;
    }
    case 'r': {
      opt.dpu_rk = atoi(optarg);
      break;
    }
    case 'l': {
      stringstream(optarg) >> opt.fld;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.sd;
      break;
    }
    case 'n': {
      opt.num = true;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'g': {
      stringstream(optarg) >> opt.gtfFile;
      break;
    }
    case 'b': {
      opt.batch_mode = true;
      opt.batch_file_name = optarg;
      break;
    }
    default: break;
    }
  }
  
  if (umi_flag) {
    opt.umi = true;
    opt.single_end = true; // UMI implies single-end reads
  }
  
  if (bus_flag) {
    if (!opt.num) {
      opt.batch_bus_write = true;
    } else {
      opt.batch_bus = true;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
 
  if (verbose_flag) {
    opt.verbose = true;
  }

  if (single_flag) {
    opt.single_end = true;
    opt.single_overhang = true;
  }
  
  if (quant_flag) {
    opt.pseudo_quant = true;
  }

  if (strand_flag) {
    opt.strand_specific = true;
  }
  
  if (strand_FR_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::FR;
  }
  
  if (strand_RF_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::RF;
  }

  if (pbam_flag) {
    opt.pseudobam = true;
  }
}

bool ParseTechnology(const std::string &techstr, BUSOptions& busopt, std::vector<std::string> &errorList) {
  auto i1 = techstr.find(':');
  if (i1 == std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), none found: \"" + techstr + "\"");    
    return false;
  }
  auto i2 = techstr.find(':', i1+1);
  if (i2 == std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), only one found: \"" + techstr + "\"");    
    return false;
  }
  auto ip = techstr.find(':', i2+1);
  if (ip != std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), three found: \"" + techstr + "\"");    
    return false;
  }
  auto bcstr = techstr.substr(0, i1);
  auto umistr = techstr.substr(i1+1,i2-i1-1);
  auto seqstr = techstr.substr(i2+1);
  


  int maxnf = 0;

  auto convert_commas_to_vector = [&](const std::string &s, std::vector<BUSOptionSubstr> &v) -> bool {
    std::vector<int> vv;
    v.clear();
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) {
      try {
        int i = stoi(t);
        vv.push_back(i);
      } catch (std::invalid_argument &e) {
        errorList.push_back("Error: converting to int: \"" + t + "\"");
        return false;
      }
    }

    int nv = vv.size();
    if (nv % 3 == 0) {
      for (int i = 0; i+2 < nv; i+=3) {
        int f = vv[i];
        int a = vv[i+1];
        int b = vv[i+2];
        if (f < 0) {
          errorList.push_back("Error: invalid file number (" + to_string(f) + ")  " + s);
        }
        if (a <  0) {
          errorList.push_back("Error: invalid start (" + to_string(a) + ")  " + s);
        }
        if (b != 0 && b <= a) {
          errorList.push_back("Error: invalid stop (" + to_string(b) + ") has to be after start (" + to_string(a) + ")  " + s);
        }
        v.push_back(BUSOptionSubstr(f,a,b));
        if (f > maxnf) {
          maxnf = f;
        }
      }
    } else {
      errorList.push_back("Error: number of values has to be multiple of 3 " + s);
      return false;
    }

    busopt.nfiles = maxnf+1;
    return true;
  };

  

  std::vector<BUSOptionSubstr> v;
  if (!convert_commas_to_vector(bcstr,v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty barcode list " + bcstr);
    return false;
  }
  busopt.bc = std::move(v);

  if (!convert_commas_to_vector(umistr, v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty UMI list " + umistr);
    return false;
  }
  busopt.umi = std::move(v);


  if (!convert_commas_to_vector(seqstr, v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty sequence list " + bcstr);
    return false;
  }

  busopt.seq = std::move(v);

  return true;
}


bool CheckOptionsPseudo(ProgramOptions& opt) {

  bool ret = true;

  cerr << endl;
  // check for index
  if (opt.index.empty()) {
    cerr << ERROR_STR << " kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  }
  
  assert(!(opt.batch_bus && opt.batch_bus_write));
  bool opt_bus_mode = opt.batch_bus || opt.batch_bus_write;

  if (opt.pseudo_quant) {
    if (!opt.batch_mode) {
      cerr << ERROR_STR << " --quant can only be run with in batch mode with --batch" << endl;
      ret = false;
    }
    if (opt_bus_mode) {
      cerr << ERROR_STR << " --quant cannot be run with --bus" << endl;
      ret = false;
    }
  }
  
  if (opt_bus_mode && opt.umi) {
    cerr << ERROR_STR << " UMI cannot be run with --bus" << endl;
    ret = false;   
  }

  // check for read files
  if (!opt.batch_mode) {
    if (opt.umi) {
      cerr << ERROR_STR << " UMI must be run in batch mode, use --batch option" << endl;
      ret = false;      
    }
    if (opt_bus_mode && ret) {
      cerr << "--bus specified; will try running read files in batch mode" << endl;
      opt.batch_ids.push_back("sample");
      opt.batch_mode = true;
      opt.pseudo_read_files_supplied = true;
      if (opt.files.size() != 1 && opt.files.size() != 2) {
        cerr << ERROR_STR << " A minimum of one and a maximum of two read files must be provided" << endl;
        ret = false;
      } else if (opt.single_end && opt.files.size() != 1) {
        cerr << ERROR_STR << " single-end mode requires exactly one read file" << endl;
        ret = false;
      } else {
        std::string f1,f2;
        struct stat stFileInfo;
        if (opt.files.size() == 1) {
          f1 = opt.files[0];
          opt.batch_files.push_back({f1});
          auto intstat = stat(f1.c_str(), &stFileInfo);
          if (intstat != 0) {
            cerr << ERROR_STR << " file not found " << f1 << endl;
            ret = false;
          }
        } else {
          f1 = opt.files[0];
          f2 = opt.files[1];
          opt.batch_files.push_back({f1,f2});
          auto intstat = stat(f1.c_str(), &stFileInfo);
          if (intstat != 0) {
            cerr << ERROR_STR << " file not found " << f1 << endl;
            ret = false;
          }
          intstat = stat(f2.c_str(), &stFileInfo);
          if (intstat != 0) {
            cerr << ERROR_STR << " file not found " << f2 << endl;
            ret = false;
          }
        }
      }
    } else if (opt.files.size() == 0) {
      cerr << ERROR_STR << " Missing read files" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;      
      for (auto& fn : opt.files) {        
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          cerr << ERROR_STR << " file not found " << fn << endl;
          ret = false;
        }
      }
    }
  } else {
    if (opt.files.size() != 0) {
      cerr << ERROR_STR << " cannot specify batch mode and supply read files" << endl;
      ret = false;
    } else {
      // check for batch files
      if (opt.batch_mode) {
        struct stat stFileInfo;
        auto intstat = stat(opt.batch_file_name.c_str(), &stFileInfo);
        if (intstat != 0) {
          cerr << ERROR_STR << " file not found " << opt.batch_file_name << endl;
          ret = false;
        }
        // open the file, parse and fill the batch_files values
        std::ifstream bfile(opt.batch_file_name);
        std::string line;
        std::string id,f1,f2;
        while (std::getline(bfile,line)) {
          if (line.size() == 0) {
            continue;
          }
          std::stringstream ss(line);
          ss >> id;
          if (id[0] == '#') {
            continue;
          }
          opt.batch_ids.push_back(id);
          if (opt.single_end && !opt.umi) {
            ss >> f1;
            opt.batch_files.push_back({f1});
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f1 << endl;
              ret = false;
            }
          } else {
            ss >> f1 >> f2;
            if (!opt.umi) {
              opt.batch_files.push_back({f1,f2});
            } else {
              opt.umi_files.push_back(f1);
              opt.batch_files.push_back({f2});
            }
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f1 << endl;
              ret = false;
            }
            intstat = stat(f2.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f2 << endl;
              ret = false;
            }
          }
        }
      }
    }
  }


  /*
  if (opt.strand_specific && !opt.single_end) {
    cerr << "Error: strand-specific mode requires single-end mode" << endl;
    ret = false;
  }*/

  if (!opt.single_end) {
    if (opt.files.size() % 2 != 0) {
      cerr << "Error: paired-end mode requires an even number of input files" << endl
           << "       (use --single for processing single-end reads)" << endl;
      ret = false;
    }
  }
  
  if (opt.umi) {
    opt.single_end = true;
    if (opt.fld != 0.0 || opt.sd != 0.0) {
      cerr << "[~warn] you supplied fragment length information for UMI data which will be ignored" << endl;
    }
  } else if (opt_bus_mode) {
    if (opt.fld != 0.0 || opt.sd != 0.0) {
      cerr << "[~warn] you supplied fragment length information for --bus mode which will be ignored" << endl;
    }
    opt.fld = 0.0;
    opt.sd = 0.0;
  } else {
    if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
      cerr << "Error: cannot supply mean/sd without supplying both -l and -s" << endl;
      ret = false;
    }

    if (opt.single_end && (opt.fld == 0.0 || opt.sd == 0.0)) {
      cerr << "Error: fragment length mean and sd must be supplied for single-end reads using -l and -s" << endl;
      ret = false;
    } else if (opt.fld == 0.0 && ret) {
      // In the future, if we have single-end data we should require this
      // argument
      cerr << "[quant] fragment length distribution will be estimated from the data" << endl;
    } else if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
      cerr << "[quant] fragment length distribution is truncated gaussian with mean = " <<
        opt.fld << ", sd = " << opt.sd << endl;
    }

    if (!opt.single_end && (opt.fld > 0.0 && opt.sd > 0.0)) {
      cerr << "[~warn] you specified using a gaussian but have paired end data" << endl;
      cerr << "[~warn] we suggest omitting these parameters and let us estimate the distribution from data" << endl;
    }
  }

  if (opt.fld < 0.0) {
    cerr << "Error: invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    cerr << "Error: invalid value for fragment length standard deviation " << opt.sd << endl;
    ret = false;
  }

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      }
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "[~warn]  you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
    if (opt.threads > 1 && opt.pseudobam) {
      cerr << "Error: pseudobam is not compatible with running on many threads."<< endl;
      ret = false;
    }
  }

  if (opt.pseudobam) {
    cerr << "Error: Pseudobam not supported yet in pseudo mode, use quant mode to obtain BAM files" << endl;
    ret = false;
  }
  
  if (!opt_bus_mode && opt.num) {
    cerr << "Warning: --bus option was not used, so --num option will be ignored" << endl;
  }
  
  if (opt_bus_mode && ret) { // Set up BUS options
    auto& busopt = opt.busOptions;
    busopt.seq.push_back(BUSOptionSubstr(0,0,0));
    busopt.umi.resize(0);
    busopt.bc.resize(0);
    if (opt.single_end) {
      busopt.nfiles = 1;
      busopt.paired = false;
    } else {
      busopt.nfiles = 2;
      busopt.seq.push_back(BUSOptionSubstr(1,0,0));
      busopt.paired = true;
    }
  }

  return ret;
}

void usage() {
  cout << "kallisto " << KALLISTO_VERSION << endl << endl
       << "Usage: kallisto <CMD> [arguments] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    index         Builds a kallisto index "<< endl
       << "    quant         Runs the quantification algorithm " << endl
       << "    quant-tcc     Runs quantification on transcript-compatibility counts" << endl
       << "    bus           Generate BUS files for single-cell data " << endl
       << "    merge         Merges several batch runs " << endl
       << "    h5dump        Converts HDF5-formatted results to plaintext" << endl
       << "    inspect       Inspects and gives information about an index" << endl 
       << "    version       Prints version information" << endl
       << "    cite          Prints citation information" << endl << endl
       << "Running kallisto <CMD> without arguments prints usage information for <CMD>"<< endl << endl;
}


void usagePseudo(bool valid_input = true) {
  if (valid_input) {
    cout << "kallisto " << KALLISTO_VERSION << endl
         << "Computes equivalence classes for reads and quantifies abundances (deprecated)" << endl << endl;
  }

  cout << "Usage: kallisto pseudo [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              pseudoalignment" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-u  --umi                     First file in pair is a UMI file" << endl
       << "-b  --batch=FILE              Process files listed in FILE" << endl
       << "    --quant                   Quantify using EM algorithm (only in batch mode)" << endl
       << "    --bus                     Output a BUS file" << endl
       << "    --single                  Quantify single-end reads" << endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << endl
       << "                              (default: -l, -s values are estimated from paired" << endl
       << "                               end data, but are required when using --single" << endl
       << "                               unless outputting a BUS file via --bus)" << endl
       << "    --fr-stranded             Strand specific reads, first read forward" << endl
       << "    --rf-stranded             Strand specific reads, first read reverse" << endl
       << "-n, --num                     Output number of read in BUS file flag column (only with --bus)" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "-d, --dpu=INT                 Number of dpus to use (default: 64)" << endl
       << "-r, --rank=INT                Number of dpu ranks to use (default: 1)" << endl;
//       << "    --pseudobam               Output pseudoalignments in SAM format to stdout" << endl;

}

std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }

  return res;
}

std::string get_local_time() {
  time_t rawtime;
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);


  if (argc < 2) {
    usage();
    exit(1);
  } else {
    auto start_time(get_local_time());
    ProgramOptions opt;
    string cmd(argv[1]);
    if (cmd == "pseudo") {
      if (argc==2) {
        usagePseudo();
        return 0;
      }
      ParseOptionsPseudo(argc-1,argv+1,opt);
      if (!CheckOptionsPseudo(opt)) {
        cerr << endl;
        usagePseudo(false);
        exit(1);
      } else {
        // pseudoalign the reads
        KmerIndex index(opt);
        index.load(opt);

        MinCollector collection(index, opt);
        int64_t num_processed = 0;
        int64_t num_pseudoaligned = 0;
        int64_t num_unique = 0;

        Transcriptome model; // empty model
        if (!opt.gtfFile.empty()) {          
          model.parseGTF(opt.gtfFile, index, opt, true);
        }
        MasterProcessor MP(index, opt, collection, model);
        std::cerr << "size of index.kmap.size_: " << index.kmap.size_ << std::endl;
        std::cerr << "size of index.kmap.size(): " << index.kmap.size() << std::endl;
        if (!opt.batch_mode) {
          num_processed = ProcessReads(MP, opt);
          collection.write((opt.output + "/pseudoalignments"));
        } else {          
          num_processed = ProcessBatchReads(MP,opt);
        }

        std::string call = argv_to_string(argc, argv);

        if (!opt.batch_bus) {
          for (int id = 0; id < MP.batchCounts.size(); id++) {
            const auto &cc = MP.batchCounts[id];
            for (const auto &p : cc) {
              if (p.first < index.num_trans) {
                num_unique += p.second;
              }
              num_pseudoaligned += p.second;
            }
          }
        } else {
          for (int i = 0; i < index.num_trans; i++) {
            num_unique += collection.counts[i];          
          }
          for (int i = 0; i < collection.counts.size(); i++) {
            num_pseudoaligned += collection.counts[i];
          }
        }
        
        std::ofstream transout_f((opt.output + "/transcripts.txt"));
        for (const auto &t : index.target_names_) {
          transout_f << t << "\n";
        }
        transout_f.close();

        plaintext_aux(
            opt.output + "/run_info.json",
            std::string(std::to_string(index.num_trans)),
            std::string(std::to_string(0)), // no bootstraps in pseudo
            std::string(std::to_string(num_processed)),
            std::string(std::to_string(num_pseudoaligned)),
            std::string(std::to_string(num_unique)),
            KALLISTO_VERSION,
            std::string(std::to_string(index.INDEX_VERSION)),
            start_time,
            call);

        
        
        std::vector<std::vector<std::pair<int32_t, double>>> Abundance_mat;
        std::vector<std::pair<double, double>> FLD_mat;
          
        if (opt.pseudo_quant) {
          int n_batch_files = opt.batch_files.size();
          Abundance_mat.resize(n_batch_files, {});
          FLD_mat.resize(n_batch_files, {});

          std::cerr << "[quant] Running EM algorithm for each cell .."; std::cerr.flush();
          
          auto EM_lambda = [&](int id) {          
            MinCollector collection(index, opt);
            collection.flens = MP.batchFlens[id];
            collection.counts.assign(index.ecmap.size(), 0);
            const auto& bc = MP.batchCounts[id];
            for (const auto &p : bc) {
              collection.counts[p.first] = p.second;
            }
            // if mean FL not provided, estimate
            std::vector<int> fld;
            if (opt.fld == 0.0) {
              fld = collection.flens; // copy
              collection.compute_mean_frag_lens_trunc(false);
            } else {
              auto mean_fl = (opt.fld > 0.0) ? opt.fld : collection.get_mean_frag_len(true);
              if (mean_fl == std::numeric_limits<double>::max()) {
                std::cerr << "Couldn't estimate fragment length for batch file number " << id << std::endl;
                return;
              }
              auto sd_fl = opt.sd;
              collection.init_mean_fl_trunc( mean_fl, sd_fl );
              fld = trunc_gaussian_counts(0, MAX_FRAG_LEN, mean_fl, sd_fl, 10000);
            }
            std::vector<int> preBias(4096,1);
            if (opt.bias) {
              //preBias = collection.bias5; // copy
            }

            auto fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);

            EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
            em.run(10000, 50, false, opt.bias);

            auto &ab_m = Abundance_mat[id];
            for (int i = 0; i < em.alpha_.size(); i++) {
              if (em.alpha_[i] > 0.0) {
                ab_m.push_back({i,em.alpha_[i]});
              }
            }

            double mean_fl = collection.get_mean_frag_len(true);
            double sd_fl = collection.get_sd_frag_len();
            FLD_mat[id] = {mean_fl, sd_fl};
          }; // end of EM_lambda

          std::vector<std::thread> workers;
          int num_ids = opt.batch_ids.size();
          int id =0;
          while (id < num_ids) {
            workers.clear();
            int nt = std::min(opt.threads, (num_ids - id));
            int first_id = id;
            for (int i = 0; i < nt; i++,id++) {
              workers.emplace_back(std::thread(EM_lambda, id));
              //workers.emplace_back(std::thread(ReadProcessor(index, opt, tc, *this, id,i)));
            }
            
            for (int i = 0; i < nt; i++) {
              workers[i].join();
            }
          }

          std::cerr << " done" << std::endl;


        }
        cerr << endl;

        std::string prefix = opt.output + "/matrix";
        std::string ecfilename = prefix + ".ec";
        std::string tccfilename = prefix + ".tcc.mtx";
        std::string abfilename = prefix + ".abundance.mtx";
        std::string cellnamesfilename = prefix + ".cells";
        std::string fldfilename = prefix + ".fld.tsv";
        std::string genelistname = prefix + ".genelist.txt";
        std::string genecountname = prefix + ".genes.mtx";
        std::string busbarcodelistname = prefix + ".barcodes";
        std::string busoutputname = opt.output + "/output.bus";

        writeECList(ecfilename, index);
        writeCellIds(cellnamesfilename, opt.batch_ids);
        if (opt.batch_bus || opt.batch_bus_write) {
          if (opt.batch_bus_write) {
            writeBUSMatrix(busoutputname, MP.batchCounts, index.ecmap.size());
          }
          if (!MP.batchCounts.empty()) {
            // Write out fake barcodes that identify each cell
            std::vector<std::string> fake_bcs;
            fake_bcs.reserve(MP.batchCounts.size());
            for (size_t j = 0; j < MP.batchCounts.size(); j++) {
              fake_bcs.push_back(binaryToString(j, BUSFORMAT_FAKE_BARCODE_LEN));
            }
            writeCellIds(busbarcodelistname, fake_bcs);
          }
          // Write out index:
          index.write((opt.output + "/index.saved"), false);
          // Write out fragment length distributions if reads paired-end:
          if (!opt.single_end) {
            std::ofstream flensout_f((opt.output + "/flens.txt"));
            for (size_t id = 0; id < opt.batch_ids.size(); id++) {
              std::vector<int> fld = MP.batchFlens[id];
              for ( size_t i = 0 ; i < fld.size(); ++i ) {
                if (i != 0) {
                  flensout_f << " ";
                }
                flensout_f << fld[i];
              }
              flensout_f << "\n";
            }
            flensout_f.close();
          }
        } else {
          writeSparseBatchMatrix(tccfilename, MP.batchCounts, index.ecmap.size());
        }
        if (opt.pseudo_quant) {
          writeSparseBatchMatrix(abfilename, Abundance_mat, index.num_trans);
          writeFLD(fldfilename, FLD_mat);
        }
        if (!opt.gtfFile.empty()) {
          // write out gene info
          std::vector<std::vector<std::pair<int32_t, double>>> geneCounts;
          geneCounts.assign(MP.batchCounts.size(), {});
          
          std::unordered_set<int> gene_ids;
          gene_ids.reserve(100);
          int n_batch_files = opt.batch_files.size();
          std::vector<double> gc;

          for (int id = 0; id < n_batch_files; id++) {
            auto& sgc = geneCounts[id];
            gc.assign(model.genes.size(), 0.0);  
            const auto& bc = MP.batchCounts[id];
            for (auto &p : bc) {
              int ec = p.first;
              if (ec < 0) {
                continue; 
              }
              if (ec < index.num_trans) {
                int g_id = model.transcripts[ec].gene_id;
                if (g_id != -1) {
                  gc[g_id] += p.second;
                }
              } else {
                gene_ids.clear();
                for (auto t : index.ecmap[ec]) {
                  int g_id = model.transcripts[t].gene_id;
                  if (g_id != -1) {
                    gene_ids.insert(g_id);
                  }
                }
                if (!gene_ids.empty()) {
                  double n_genes = gene_ids.size();
                  for (auto &g_id : gene_ids) {
                    gc[g_id] += p.second / n_genes;
                  }
                }
              }
            }

            for (int j = 0; j < gc.size(); j++) {
              if (gc[j] > 0.0) {
                sgc.push_back({j, gc[j]});                
              }
            }
          }


         

          writeGeneList(genelistname, model);
          writeSparseBatchMatrix(genecountname, geneCounts, model.genes.size());
        }


        if (opt.pseudobam) {       
          std::vector<double> fl_means(index.target_lens_.size(),0.0);
          EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
          MP.processAln(em, false);
        }
      }

      
    }  
    else {
      cerr << "Error: invalid command " << cmd << endl;
      usage();
      exit(1);
    }

  }

  fflush(stdout);

  return 0;
}
