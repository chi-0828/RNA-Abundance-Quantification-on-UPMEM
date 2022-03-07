#ifndef DPU_DEF_H
#define DPU_DEF_H

#define MAX_READ_LEN 150
#define MAX_READ_N 800
#define MAX_PACKET_SIZE (MAX_READ_LEN*MAX_READ_N)

#define MAX_table_n 4000000

typedef struct dpu_args {
    // k-mer size
    int32_t k;
    // mapping or intersection
    int32_t run;
}dpu_args;

#endif