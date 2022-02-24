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



