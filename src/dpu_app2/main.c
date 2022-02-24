#include <mram.h>
#include <assert.h>
#include "defs.h"
#include "mutex.h"
#include <string.h>

#define MAX_entry 4000
#define MAX_READ_N 25
#define MAX_vector_len 4000

// non change DB
__mram_noinit int32_t vector[MAX_entry*MAX_vector_len];
__mram_noinit int32_t vector_size[MAX_entry];
__mram_noinit int32_t vector_head[MAX_entry];

// input
__host int32_t matched_id[MAX_entry];
__host int32_t matched_count[MAX_READ_N];
__host int32_t read_n;

// output
__mram_noinit int32_t result[MAX_READ_N*MAX_vector_len];

int roundTo(int value)
{
    return (value + (8 - 1)) & ~(8 - 1);
}

int main(){

  int matched_id_head = 0;
  int u[MAX_vector_len];
  int u_size = 0;
  for(int read_id = 0; read_id < read_n; read_id++){
    for(int count = 0 ; count < matched_count[read_id]; count++){
      int id = matched_id[count + matched_id_head];
      int head = vector_head[id];
      int size = vector_size[id];
      __dma_aligned int vector_cache[MAX_vector_len];
      mram_read(vector+head, vector_cache, roundTo(size));

      // intersection
      int u_tmp[MAX_vector_len];
      int u_tmp_len = 0;
      for(int a = 0, b = 0; a < u_size && b < size; ){
        if (u[a] < vector_cache[b]) {
          ++a;
        } else if (vector_cache[b] < u[a]) {
          ++b;
        } else {
          // match
          u_tmp[u_tmp_len++] = u[a];
          ++a;
          ++b;
        }
      }
      u_size = u_tmp_len;
      memcpy(u,u_tmp,u_size);
    }
    matched_id_head += matched_count[read_id];
  }
  
}