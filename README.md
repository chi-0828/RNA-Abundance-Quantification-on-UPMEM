# RNA-Abundance-Quantification-on-UPMEM

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
E.g.,
``` shell
time ./D_kallisto pseudo -i ~/data/experiment/11-mer.idx -o out --single ~/data/experiment/RNA_read/100K.fastq -l 150 -s 30 -t 8 -d 64
```

## more information
DPU code is in src/dpu_app
host code is mainly in src/ProcessReads.cpp

## Reference
### kallisto
https://github.com/pachterlab/kallisto
### UPMEM
https://github.com/CMU-SAFARI/prim-benchmarks<br>
https://sdk.upmem.com/2021.3.0/

### Testing
time ./kallisto pseudo -i ~/data/experiment/7-mer.idx -o out --single ~/data/experiment/RNA_read/100K.fastq -l 150 -s 30 -t 1



