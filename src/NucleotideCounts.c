#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "SynExtend.h"
#include "SEutils.h"

int decodeChar(char c){
  int val = -1;
  switch(c){
    case 'C':
    case 'c':
      val++;
    case 'G':
    case 'g':
      val++;
    case 'T':
    case 't':
      val++;
    case 'A':
    case 'a':
      val++;
    default:
      return val;
  }
}

int oligotoIndex(char *nuc, int length){
  int out = 0;
  int v;
  for (int i=0; i<length; i++){
    v = pow(4,i);
    v = v*decodeChar(*(nuc+i));
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

char *removeGaps(const unsigned char *nuc, uint64_t len){
  // Needs to be one larger since R strings aren't null-terminated
  char *newstring = calloc(len+1, sizeof(char));
  char *slow = newstring;
  const unsigned char *fast = nuc;
  fast = nuc;
  for (uint64_t i=0; i<len; i++){
    switch(*(fast+i)){
      case 'C':
      case 'c':
      case 'G':
      case 'g':
      case 'T':
      case 't':
      case 'A':
      case 'a':
      default: //comment out for no remove gaps
        *(slow++) = (char) *(fast+i);
    }
  }
  *slow = '\0';
  //Rprintf("%lu %lu\n", len+1, slow-newstring+1);
  return(newstring);
}

SEXP StringToNVDC(SEXP DNASTRING, SEXP REMOVEGAPS, SEXP EXTENDED){
  uint64_t sl = XLENGTH(DNASTRING);
  bool extended = LOGICAL(EXTENDED)[0];
  bool rgaps = LOGICAL(REMOVEGAPS)[0];
  int length = extended ? 92 : 12;
  int iterlen = extended ? 3 : 1;
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
    dnastring = removeGaps(instring, sl);  
  } else {
    dnastring = calloc(sl+1, sizeof(char));
    for(uint64_t i=0; i<sl; i++)
      dnastring[i] = (char) instring[i];
    dnastring[sl] = '\0';
  }
  

  int idx;
  for(uint64_t i=0; i<sl; i++){
    for (int j=0; j<iterlen; j++){
      if (i+j<sl){
        idx = oligotoIndex(dnastring+i, j+1);
        // increment raw count
        if (idx >= 0)
          outvals[idx] += 1.0;
        // increment mean distance
        if (j==0)
          outvals[idx+4] += (i+1.0);
      } 
      else break;
    }
  }

  // Fix means
  for (int i=0;i<4;i++)
    if (outvals[i] != 0)
      outvals[i+4] /= outvals[i];

  // Get normalized second moment
  for (uint64_t i=0; *(dnastring+i); i++){
    idx = decodeChar(*(dnastring+i));
    if (idx >= 0)
      outvals[idx+8] += pow((i+1) - outvals[idx+4], 2) / (outvals[idx]*sl);
  }

  SEXP retval = PROTECT(allocVector(REALSXP, length));
  memcpy(REAL(retval), outvals, length*sizeof(double));
  free(outvals);
  free(dnastring);
  UNPROTECT(1);
  return(retval);
}