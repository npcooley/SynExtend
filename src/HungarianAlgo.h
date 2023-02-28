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

// Global variables for freeing
bool has_alloced_mem, has_alloced_vec, has_alloced_assign;
uint8_t *sa;
bool *cc, *cr;
int *av;
double *vec;

// Steps
int* hungarian(double *costMatrix, int n);
double* correctMat(double *mat, int nrow, int ncol);
void hg_step1(double *costMatrix, int n);
void hg_step2(double *costMatrix, int n, int* assignment);
void hg_step3(double *costMatrix, uint8_t *statusarray, int n);
int hg_step4(uint8_t *statusarray, bool *coveredCol, bool *coveredRow, int n, int statuscode);
void hg_step5(double *costMatrix, bool *coveredCol, bool *coveredRow, int n);
bool hg_statuscheck(uint8_t *statusarray, int *assignment, int n);