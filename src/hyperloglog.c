#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <R_ext/Altrep.h>
#include <R_ext/Itermacros.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

#include "SEutils.h"

/*
 * HyperLogLog estimation of cardinality of large datafiles
 * Needs a hash function defined for any arbitrary input.
 * I'm probably just going to start with the simple ones,
 * INTSXP, CHARSXP, LGLSXP, STRSXP, RAWSXP
 * Anything else will just throw a not implemented error.
 *
 * Much of the random number generation was adapted from 
 * https://prng.di.unimi.it/xoshiro256starstar.c
 */

/* Hash Functions */
static uint64_t hashi64(uint64_t h) {
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdL;
  h ^= h >> 33;
  h *= 0xc4ceb9fe1a85ec53L;
  h ^= h >> 33;
  return h;
}

static inline uint64_t firstbbits(uint64_t i, int b){
  return i >> (64-b);
}

static inline int log2i(uint64_t i, int maxv){
  int j=-1;
  while(i >>= 1 && j<maxv) j++;
  return j;
}

static inline int max(int a, int b){
  return a > b ? a : b;
}

static inline int count_zeros(uint64_t i){
  int j=0;
  while(!(i&1)){
    j+=1;
    i >>=1;
  }

  return j;
}

static uint64_t HLLint(int* v, int length, int m, int nruns, int nthreads){
  // m must be a power of two, if not we'll round up
  if ((m&(m-1) || !m)) 
    m = 1 << (log2i(m,64)+1);

  int registers[m];
  for(int i=0; i<m; i++)
    registers[i] = 0;
  int h;
  // iterate over all the elements, storing max number of ending zeros in each register
  // registers are indexed by the first b bits
  #pragma omp parallel for private(h) shared(registers) num_threads(nthreads)
  for(int i=0; i<length; i++){
    h = hashi64(v[i]);
    registers[h & !(m-1)] = max(registers[h & !(m-1)], count_zeros((h & m-1) | m));
  }

  // now registers has the max counts for each bin
  float harm_mean = 0;
  int num_zero = 0;
  for(int i=0; i<m; i++){
    harm_mean += 1 / (1 << registers[i]);
    num_zero += registers[i] == 0;
  }

  // linear counting for small cardinalities
  if(num_zero > 0)
    return (int) m * log2i(m/num_zero, 64);

  harm_mean = 1 / harm_mean;
  float aval = 0;
  switch(m){
    case 0:
    case 1:
    case 2:
    case 4: 
    case 8:
    case 16:
      aval = 0.673;
      break;
    case 32:
      aval = 0.697;
      break;
    case 64:
      aval = 0.709;
      break;
    default:
      aval = 0.7213 / (1+1.079/m);
  }

  return (int) aval * m * m * harm_mean;
}