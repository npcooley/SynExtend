#include "SEutils.h"

int *sample(int n){
  int *r = malloc(sizeof(int) * n);
  int j;
  for (int i=0; i<n; i++){
    j = irand() % (i+1);
    if (j != i)
      r[i] = r[j];
    r[j] = i;
  }

  return r;
}

void shuffle_int_(int *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void shuffle_uint_(uint *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}


void shuffle_double_(double *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void shuffle_char_(char *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

/* needs C11 :/
#define shuffle(x, n) _Generic((x),
  int *: shuffle_i(x, n),
  double *: shuffle_d(x, n),
  char *: shuffle_c(x, n))
*/