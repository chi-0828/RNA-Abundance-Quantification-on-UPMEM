# RNA-Abundance-Quantification-on-UPMEM
- More details can be found in the paper: "RNA-seq Quantification on Processing in Memory Architecture: Observation and Characterization"
- [paper](https://doi.org/10.1109/NVMSA56066.2022.00014)
  
## Our new project about RNA-seq Quantification on UPMEM DPU 
[UpPipe](https://github.com/chi-0828/UpPipe)

## Cite the paper if you use D_kallisto in your work
> Liang-Chi Chen, Shu-Qi Yu, Chien-Chung Ho, Yuan-Hao Chang, Da-Wei Chang, Wei-Chen Wang, Yu-Ming Chang,
> "RNA-seq Quantification on Processing in memory Architecture: Observation and Characterization,"
> The 11th IEEE Non-Volatile Memory Systems and Applications Symposium (NVMSA), August 23-25, 2022
```
@inproceedings{chen2022rna,
  title={RNA-seq Quantification on Processing in memory Architecture: Observation and Characterization},
  author={Chen, Liang-Chi and Yu, Shu-Qi and Ho, Chien-Chung and Chang, Yuan-Hao and Chang, Da-Wei and Wang, Wei-Chen and Chang, Yu-Ming},
  booktitle={2022 IEEE 11th Non-Volatile Memory Systems and Applications Symposium (NVMSA)},
  pages={26--32},
  year={2022},
  organization={IEEE}
}
```


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
DPU allocation and CPU-DPU(DPU-CPU) transfers are in src/ProcessReads.cpp<br>

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



