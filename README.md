# RNA-Abundance-Quantification-on-UPMEM
## More destails can be found in the paper
https://scholar.google.com.tw/citations?view_op=view_citation&hl=zh-TW&user=SoyMWUsAAAAJ&citation_for_view=SoyMWUsAAAAJ:u-x6o8ySG0sC

## Build
``` shell
// build htslib first
cd ext/htslib
autoheader
autoconf
make -j16
// build our main program
cd ../..
mkdir obj
cd src 
make -j16
```

## Usage
``` shell
./D_kallisto pseudo [fastq file] 
      -i [index file] 
      -o [output path] 
      -t [num of CPU threads] 
      -d [num DPUs]
      --single
      -l [double]
      -s [double]
```
E.g., testing 100K reads/11-mer by 64*8 dpus 
``` shell
time ./D_kallisto pseudo -i ~/data/experiment/11-mer.idx -o out --single ~/data/experiment/RNA_read/100K.fastq -l 150 -s 30 -t 8 -d 64
```

## More information
DPU program is in src/dpu_app<br>
DPU allocation and CPU-DPU(DPU-CPU) transfers is in src/ProcessReads.cpp<br>

## Reference
### kallisto
https://github.com/pachterlab/kallisto
### UPMEM
https://github.com/CMU-SAFARI/prim-benchmarks<br>
https://sdk.upmem.com/2021.3.0/<br>
https://sdk.upmem.com/2021.3.0/CppAPI/index.html<br>

### Testing CPU-based kallisto
``` shell
time ./kallisto pseudo -i ~/data/experiment/11-mer.idx -o out --single ~/data/experiment/RNA_read/100K.fastq -l 150 -s 30 -t 8 
```



