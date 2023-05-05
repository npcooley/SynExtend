#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

#include "SynExtend.h"
#include "SEutils.h"

#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)

SEXP pseudoRandomSample(SEXP N){
  const uint64_t seed = 0xAE356E5F366847A2; //INTEGER(SEED)[0];
  int n = INTEGER(N)[0];
  

  struct RNGstate64 *r = malloc(sizeof(struct RNGstate64));
  seedRNGState64(r, seed);


  SEXP outvec = PROTECT(allocVector(INTSXP, n));
  int *ptr = INTEGER(outvec);
  for (int i=0; i<n; i++)
    ptr[i] = (xorshift128p(r) >> 32);

  free(r);
  UNPROTECT(1);
  return outvec;
}

SEXP seededPseudoRandomSample(SEXP N, SEXP SEED){
  if (LENGTH(SEED) < 2){
    error("SEED must be an integer vector of length 2\n");
  }
  const uint64_t seed = (((uint64_t) INTEGER(SEED)[0]) << 32) | INTEGER(SEED)[1];
  int n = INTEGER(N)[0];
  
  struct RNGstate64 *r = malloc(sizeof(struct RNGstate64));
  seedRNGState64(r, seed);

  SEXP outvec = PROTECT(allocVector(INTSXP, n));
  int *ptr = INTEGER(outvec);
  for (int i=0; i<n; i++)
    ptr[i] = (xorshift128p(r) >> 32);

  free(r);
  UNPROTECT(1);
  return outvec;
}

SEXP randomProjection(SEXP VEC, SEXP NONZERO, SEXP N, SEXP OUTDIM, SEXP NTHREADS){
  int64_t m = INTEGER(OUTDIM)[0];
  double *v = REAL(VEC);
  int threads = INTEGER(NTHREADS)[0];
  const int bitwidth = 64;
  const int remainder = m % bitwidth;

  int num_nonzero = INTEGER(N)[0];
  int *nzpos = INTEGER(NONZERO);

  struct RNGstate64 *r = malloc(sizeof(struct RNGstate64));
  uint64_t randnumber;
  int idx, loc;

  // new distribution: +/- 1 with probability 1/2
  SEXP outvec = PROTECT(allocVector(REALSXP, m));
  double *outvals = REAL(outvec);
  memset(outvals, 0, m*sizeof(double));

  double val, invval;
  int i,j,k;
  #pragma omp parallel for private(i, j, k, val, invval) num_threads(threads)
  for (i=0; i<num_nonzero; i++){
    loc = nzpos[i]-1;
    val = v[loc];
    invval = -1 * val;
    seedRNGState64(r, (uint64_t) loc);
    idx = 0;
    for (j=0; j < (m/bitwidth); j++){
      randnumber = xorshift128p(r);
      for (k=0; k<bitwidth; k++)
        outvals[idx+k] += randnumber & (1ULL << k) ? val : invval; 
      idx += bitwidth;
    }

    // faster to check once than to roll another number
    randnumber = remainder==0 ? 0 : xorshift128p(r);
    for (k=0; k<remainder; k++)
        outvals[idx+k] += randnumber & (1ULL << k) ? val : invval; 
  }

  free(r);
  UNPROTECT(1);
  return outvec;
}

// in-place addition of cophenetic distances
// This is an exact replica of DECIPHER code, just including for namespacing and ease
// code is short enough to include explicitly
SEXP se_cophenetic(SEXP Index1, SEXP Index2, SEXP N, SEXP D, SEXP H)
{
  int i, j, val;
  int *I = INTEGER(Index1);
  int *J = INTEGER(Index2);
  int l1 = length(Index1);
  int l2 = length(Index2);
  int n = asInteger(N);
  double *d = REAL(D);
  double h = asReal(H);
  
  for (i = 0; i < l1; i++) {
    for (j = 0; j < l2; j++) {
      if (I[i] < J[j]) {
        val = n*(I[i] - 1) - I[i]*(I[i] - 1)/2 + J[j] - I[i] - 1;
      } else {
        val = n*(J[j] - 1) - J[j]*(J[j] - 1)/2 + I[i] - J[j] - 1;
      }
      d[val] += h;
    }
  }
  
  return D;
}