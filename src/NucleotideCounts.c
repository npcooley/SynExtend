#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// benchmarking
#include <time.h>


// for OpenMP parallel processing
#ifdef _OPENMP
// see https://github.com/Bioconductor/SparseArray/issues/9#issuecomment-2096839890
#undef match
#include <omp.h>
#endif


#include "SynExtend.h"
#include "SEutils.h"

clock_t startt, endt;
double cpu_time_used;

int decodeChar(char c, const char *lookup){
  // convert to lowercase
  char newc = c > 96 ? c-32 : c;
  for(int i=0;lookup[i];i++)
    if (lookup[i] == newc)
      return i;
  return -1;
}

int oligotoIndex(char *nuc, int length, const char *codex, int clen){
  int out = 0;
  int v;
  for (int i=0; i<length; i++){
    v = pow(clen,i);
    v = v*decodeChar(*(nuc+i), codex);
    if (v < 0)
      return -1;
    out += v;
  }
  if (length==2)
    out += 12;
  else if (length==3)
    out += 28;
  return out;
}

char *removeGaps(const unsigned char *nuc, uint64_t len, const char *codex){
  // Needs to be one larger since R strings aren't null-terminated
  char *newstring = calloc(len+1, sizeof(char));
  char *slow = newstring;
  const unsigned char *fast = nuc;
  fast = nuc;
  for (uint64_t i=0; i<len; i++)
    if (decodeChar(fast[i], codex) >= 0)
      *(slow++) = (char) fast[i];
  *slow = '\0';

  return(newstring);
}

SEXP StringToNVDT(SEXP DNASTRING, SEXP REMOVEGAPS, SEXP EXTENDED, SEXP USEDNA){
  uint64_t sl = XLENGTH(DNASTRING);
  bool extended = LOGICAL(EXTENDED)[0];
  bool rgaps = LOGICAL(REMOVEGAPS)[0];

  // DNA or AA setup
  bool useDNA = LOGICAL(USEDNA)[0];
  int length, codexlen;
  const char *codex;
  if (useDNA){
    length = extended ? 92 : 12;
    codexlen = 4;
    codex = "ATGC";
  } else {
    length = 60;
    codexlen = 20;
    codex = "ARNDCQEGHILKMFPSTWYV";
  }

  int iterlen = (useDNA & extended) ? 3 : 1;
  /*
   * outvals is a 12 or 92 length vector
   * raw counts stored as follows:
   *
   *   A   T   G   C
   *  AA  AT  AG  AC
   *  TA  TT  TG  TC
   *       ...
   * AAA AAT AAG AAC
   * ATA ATT ATG ATC
   *       ...
   *
   * Indexing: A=0, T=1, G=2, C=3
   * Dinucleotides: 12 + 4*N2 + N1
   * Trinucleotides: 28 + 16*N3 + 4*N2 + N1
   */

  // Get Raw Counts
  double *outvals = calloc(length, sizeof(double));
  const unsigned char *instring = RAW(DNASTRING);
  char *dnastring;
  if (rgaps){
    dnastring = removeGaps(instring, sl, codex);
  } else {
    dnastring = calloc(sl+1, sizeof(char));
    for(uint64_t i=0; i<sl; i++)
      dnastring[i] = (char) instring[i];
    dnastring[sl] = '\0';
  }


  int idx;
  uint64_t i, j;
  for(i=0; i<sl; i++){
    for (j=0; j<iterlen; j++){
      if (i+j<sl){
        idx = oligotoIndex(dnastring+i, j+1, codex, codexlen);
        // increment raw count
        if (idx >= 0){
          outvals[idx] += 1.0;
          // increment mean distance
          if (j==0)
            outvals[idx+codexlen] += (i+1.0);
        }
      }
      else break;
    }
  }

  // Fix means
  for (int i=0;i<codexlen;i++)
    if (outvals[i] != 0)
      outvals[i+codexlen] /= outvals[i];

  // Get normalized second moment
  for (uint64_t i=0; *(dnastring+i); i++){
    idx = decodeChar(*(dnastring+i), codex);
    if (idx >= 0)
      outvals[idx+2*codexlen] += pow((i+1) - outvals[idx+codexlen], 2) / (outvals[idx]*sl);
  }

  SEXP retval = PROTECT(allocVector(REALSXP, length));
  memcpy(REAL(retval), outvals, length*sizeof(double));
  free(outvals);
  free(dnastring);
  UNPROTECT(1);
  return(retval);
}

SEXP fastPearsonC(SEXP V1, SEXP V2){
  // assume we've already trimmed NAs
  int l1 = LENGTH(V1);
  int l2 = LENGTH(V2);
  // going to the first value, ideally they'll be the same length
  int l = l1 <= l2 ? l1 : l2;
  double *v1 = REAL(V1);
  double *v2 = REAL(V2);

  double sumx=0, sumy=0, sumprod=0;
  double sumxsq=0, sumysq=0;
  double x,y,n=0;

  for(int i=0; i<l; i++){
    x = v1[i];
    y = v2[i];
    if(ISNA(x) || ISNA(y)) continue;
    sumx += x;
    sumy += y;
    sumprod += x*y;
    sumxsq += x*x;
    sumysq += y*y;
    n++;
  }
  double r,t;

  if (n != 0){
    // calculate Pearson's correlation
    r = (n*sumprod - sumx*sumy) / (sqrt(n*sumxsq - sumx*sumx) * sqrt(n*sumysq - sumy*sumy));
    t = r * sqrt((n-2) / (1-r*r));
  } else {
    r = 0;
    t = 0;
  }



  SEXP outval = PROTECT(allocVector(REALSXP, 3));
  REAL(outval)[0] = r;
  REAL(outval)[1] = t;
  REAL(outval)[2] = n;

  UNPROTECT(1);
  return(outval);
}

static double logBaseVal;
static inline double fastlog2(double x){ return( x == 0 ? 0 : log2(x)); }

static double calc_mi_colpair(int* col1, int* col2, int collen, double addVal, int u1, int u2, double base){
  double* contTable = calloc(u1*u2, sizeof(double));
  double jEnt = 0;
  //double pseudocount = 0.05;
  double total=0;

  for(int i=0; i<collen; i++){
      contTable[col1[i]*u2+col2[i]]+=1;
      total += 1;
  }

  for(int i=0; i<u1*u2;i++){
    contTable[i] = (contTable[i]) / total;
    jEnt -= contTable[i] * (fastlog2(contTable[i])) * logBaseVal;
  }

  free(contTable);
  return(jEnt);
}

static double calc_ent_col(int* col, int collen, double addVal, int llen, double base){
  double* ptab = calloc(llen, sizeof(double));
  double ent = 0;
  //double pseudocount = 0.05;
  double total=0;
  for(int i=0; i<collen; i++){
      ptab[col[i]]+=1;
      total++;
  }

  for(int i=0; i<llen; i++){
    // could add pseudocounts here
    ptab[i] = (ptab[i]) / total;
    ent -= ptab[i] * fastlog2(ptab[i]) * logBaseVal;
  }

  free(ptab);
  return(ent);
}

SEXP MIForSequenceSets(SEXP M1, SEXP M2, SEXP NSEQS, SEXP U1, SEXP U2, SEXP BASE, SEXP NTHREADS){
  /*
   * inputs:
   *  - two matrices of residues/bases such encoded as integers
   *  - number of sequences
   *  - number of columns for each matrix (is this necessary?)
   *  - base of log used for calculations
   *
   * Note that we add an extra space for gaps, but we don't utilize it
   * This ensures we count enough values, but we don't need to
   */

  double logBase = REAL(BASE)[0];
  int nseqs = INTEGER(NSEQS)[0];
  int *m1 = INTEGER(M1);
  int *m2 = INTEGER(M2);
  int c1 = LENGTH(M1) / nseqs;
  int c2 = LENGTH(M2) / nseqs;
  int u1 = INTEGER(U1)[0];
  int u2 = INTEGER(U2)[0];
  int nthreads = INTEGER(NTHREADS)[0];
  double addVal = 1.0/((double)nseqs);

  double *ent1 = calloc(c1, sizeof(double));
  double *ent2 = calloc(c2, sizeof(double));
  double *jointEnt = calloc(c1*c2, sizeof(double));

  logBaseVal = 1.0/log2(logBase);

  for(int i=0;i<c1;i++)
    ent1[i] = calc_ent_col(m1+nseqs*i, nseqs, addVal, u1, logBase);

  for(int i=0;i<c2;i++)
    ent2[i] = calc_ent_col(m2+nseqs*i, nseqs, addVal, u2, logBase);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads)
  #endif
  for(int i=0; i<c1; i++)
    for(int j=0; j<c2; j++)
      jointEnt[i*c2+j] = calc_mi_colpair(m1+nseqs*i, m2+nseqs*j, nseqs, addVal, u1, u2, logBase);

  for(int i=0; i<c1; i++)
    for(int j=0; j<c2; j++)
      //jointEnt[i*c2+j] = jointEnt[i*c2+j] ? ((ent1[i] + ent2[j] - jointEnt[i*c2+j]) / jointEnt[i*c2+j]) : 0;
      jointEnt[i*c2+j] = ent1[i] + ent2[j] - jointEnt[i*c2+j];

  SEXP rval = PROTECT(allocVector(REALSXP, c1*c2));
  memcpy(REAL(rval), jointEnt, sizeof(double)*c1*c2);

  free(ent1);
  free(ent2);
  free(jointEnt);

  UNPROTECT(1);
  return(rval);
}
