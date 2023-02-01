#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "SynExtend.h"

double calcMoranVar(double **weights, double weightsum, double *vals, double sampmean, int n, double expval){
  // lots of values to calculate here
  // see https://en.wikipedia.org/wiki/Moran%27s_I#Expected_value
  // I'm rolling as many calculations as I can into a single loop, since it's iteration-heavy
  double s1,s2,s3,s4,s5;

  // S1, S2, S3
  s1 = 0;
  s2 = 0;
  double tmp1, tmp2, tmp3, numer=0, denom=0;
  for (int i=0; i<n; i++){
    tmp2 = 0;
    for (int j=0; j<n; j++){
      tmp1 = weights[i][j] + weights[j][i];
      tmp2 += tmp1;
      s1 += tmp1*tmp1;
    }
    s2 += tmp2*tmp2;
    tmp3 = vals[i] - sampmean;
    numer += pow(tmp3, 4);
    denom += pow(tmp3, 2);
  }
  s1 /= 2;

  // Cleaning up S3
  numer /= n;
  denom = pow((denom / n), 2);

  if (denom == 0){
    return 0;
  }
  s3 = numer / denom;

  // S4, S5
  double w2 = weightsum * weightsum;
  double n2 = n * n;
  s4 = (n2 - 3*n + 3) * s1 - n*s2 + 3*w2;
  s5 = (n2 - n)*s1 - 2*n*s2 + 6*w2;

  numer = n*s4 - s3*s5;
  denom = (n-1)*(n-2)*(n-3)*w2;

  if(denom == 0){
    return 0;
  }

  return(numer/denom - expval*expval);
}

SEXP MoransI(SEXP VALS, SEXP DIST, SEXP DIM){
  // Read in values from R
  double *vals = REAL(VALS);
  double *dist = REAL(DIST);
  int dim = INTEGER(DIM)[0];

  // Allocate some space, convert Dist object to a matrix
  double *rowsums = calloc(dim, sizeof(double));
  double **newdist = calloc(dim, sizeof(double*));
  for(int i=0; i<dim; i++) 
    newdist[i] = calloc(dim, sizeof(double));
  
  int ctr=0;
  double cur;
  for (int i=0; i<dim; i++){
    for (int j=i+1; j<dim; j++){
      cur = dist[ctr];
      newdist[i][j] += cur;
      newdist[j][i] += cur;
      rowsums[i] += cur;
      rowsums[j] += cur;
      ctr++;
    }
  }

  // Row normalize and get sample mean
  double weightsum = 0;
  double sampmean = 0;
  for (int i=0; i<dim; i++){
    for (int j=0; j<dim; j++){
      newdist[i][j] /= rowsums[i];
      weightsum += newdist[i][j];
    }
    sampmean += vals[i];
  }
  sampmean /= dim;
  
  // Calculate denominator
  double denominator = 0, numerator = 0;
  double tmp;
  for (int i=0; i<dim; i++){
    tmp = vals[i] - sampmean;
    denominator += tmp*tmp;
  }
  
  SEXP retval = R_NilValue;
  bool worked = false;
  if (denominator != 0 && weightsum!=0){
    double vi, vj;
    // calculate numerator
    for (int i=0; i<dim; i++){
      vi = vals[i];
      for (int j=0; j<dim; j++){
        vj = vals[j];
        numerator += newdist[i][j] * (vi - sampmean) * (vj - sampmean); 
      }
    }
    double moran, expected, variance;
    moran = (((double)dim) / weightsum) * (numerator / denominator);
    expected = -1.0 / (dim - 1);
    variance = calcMoranVar(newdist, weightsum, vals, sampmean, dim, expected);
    retval = PROTECT(allocVector(REALSXP, 3));
    REAL(retval)[0] = moran;
    REAL(retval)[1] = expected;
    REAL(retval)[2] = variance;
    worked = true;
  }
  
  // freeing values
  for (int i=0; i<dim; i++){
    free(newdist[i]);
  }
  free(newdist);
  free(rowsums);
  if (worked)
    UNPROTECT(1);
  return retval;
}








