#ifndef FORTRAN_INCLUDE
#define FORTRAN_INCLUDE

#include <R.h>

void F77_SUB(GetRNGstate)(void) { GetRNGstate(); }
void F77_SUB(PutRNGstate)(void) { PutRNGstate(); }
double F77_SUB(unif_rand)(void) { return unif_rand(); }
double F77_SUB(norm_rand)(void) { return norm_rand(); }

extern void F77_NAME(tabulate_double)(double *v, int *l, double *out_val, int *out_count, int *ctr);
extern void F77_NAME(tabulate_int)(int *v, int *l, int *out_val, int *out_count, int *ctr);
extern void F77_NAME(find_gini_split)(double *v, int *response, int *l, int *nclass, double *o_v, double *o_gini_score);
extern void F77_NAME(find_sse_split)(double *v, double *response, int *l, double *o_v, double *o_sse_score);

#endif