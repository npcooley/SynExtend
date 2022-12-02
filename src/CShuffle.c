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

/*** R .C Functions ***/
void shuffleRInt(int *v, int *l){ 
  GetRNGstate();
  shuffle(int, v, *l); 
  PutRNGstate();
}

void shuffleRRepl(int *v, int *l){ 
  GetRNGstate();
  struct RNGstate32 *r = malloc(sizeof(struct RNGstate32));
  uint64_t seed = (((uint64_t)irand()) << 32) | ((uint64_t)irand());
  seedRNGState32(r, seed);
  for (int i=0; i<*l; i++)
    v[i] = (xorshift32b(r) >> 1);

  PutRNGstate();
  return;
}


