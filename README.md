# RNA-Abundance-Quantification-on-UPMEM

## build
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

### DPU code is in src/dpu_app
### host code is mainly in src/ProcessReads.cpp

### Reference
#### https://github.com/pachterlab/kallisto

### Testing
time ./kallisto pseudo -i ~/data/experiment/7-mer.idx -o out --single ~/data/experiment/RNA_read/100K.fastq -l 150 -s 30 -t 1



