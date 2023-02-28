#include "HungarianAlgo.h"

SEXP HungarianAssignment(SEXP MATVEC, SEXP DIM){
  /*
   * Arguments:
   * - MATVEC: numeric vector containing costs. Should be non-negative
   * -    DIM: integer vector containing (nrow, ncol)
   */
  has_alloced_mem = false;
  has_alloced_vec = false;
  has_alloced_assign = false;
  sa = NULL;
  cc = NULL;
  cr = NULL;
  int nrow = INTEGER(DIM)[0];
  int ncol = INTEGER(DIM)[1];

  // correctMat copies the object, no inplace modification
  vec = correctMat(REAL(MATVEC), nrow, ncol);
  has_alloced_vec = true;

  int maxdim = nrow > ncol ? nrow : ncol;

  // res is allocated in hungarian()
  int *res = hungarian(vec, maxdim);

  SEXP outval = PROTECT(allocVector(INTSXP, maxdim));
  memcpy(INTEGER(outval), res, maxdim*sizeof(int));

  free(res);
  has_alloced_assign = false;
  free(vec);
  has_alloced_vec = false;

  UNPROTECT(1);
  return(outval);
}

void hungarianCleanup(){
  if(has_alloced_mem){
    free(cc);
    free(cr);
    free(sa);
  }

  if(has_alloced_assign){
    free(av);
  }

  if(has_alloced_vec){
    free(vec);
  }
}

int *hungarian(double *costMatrix, int n){
  hg_step1(costMatrix, n);

  av = calloc(n, sizeof(int));
  has_alloced_assign = true;

  hg_step2(costMatrix, n, av);

  return av;
}


double *correctMat(double *mat, int nrow, int ncol){
  double *newmat = NULL;
  if (nrow > ncol) {
    // in this case, we have to add more columns then zero them
    newmat = (double *)calloc(nrow*nrow, sizeof(double));
    int ctr = nrow*ncol-1;
    for(int i=nrow-1; i>=0; i--){
      for (int j=nrow-1;j>=0; j--){
        if (j<ncol) newmat[i*nrow+j] = mat[ctr--];
      }
    }
  } else {
    // have to add rows then zero them out
    newmat = (double *)calloc(ncol*ncol, sizeof(double));
    memcpy(newmat, mat, nrow*ncol*sizeof(double));
  }

  return newmat;
}

void hg_step1(double *costMatrix, int n){
  // Step 1: subtract minimum of each row from each row,
  //         and repeat for each column

  // Rows
  double minv = -1;
  int idx;
  for (int i=0; i<n; i++){
    idx = n*i;
    minv = costMatrix[idx];
    for (int j=1; j<n; j++){
      if (minv == 0) break;
      minv = minv > costMatrix[idx+j] ? costMatrix[idx+j] : minv; 
    }

    if (minv != 0)
      for (int j=0; j<n; j++)
        costMatrix[idx+j] -= minv;
  }

  // Columns
  for (int i=0; i<n; i++){
    idx = i;
    minv = costMatrix[idx];
    for (int j=1; j<n; j++){
      if (minv == 0) break;
      minv = minv > costMatrix[idx+j*n] ? costMatrix[idx+j*n] : minv; 
    }

    if (minv != 0)
      for (int j=0; j<n; j++)
        costMatrix[idx+j*n] -= minv;
  }

  return;
}

void hg_step2(double *costMatrix, int n, int* assignment){
  // status array holds the status of each position as 4 bit string
  // abcd, a=starred, b=primed, c=zero, d=covered 
  //uint8_t *statusarray = calloc(n*n, 1);
  //bool *coveredCol, *coveredRow;
  sa = calloc(n*n, 1);
  cc = calloc(n, sizeof(bool));
  cr = calloc(n, sizeof(bool));
  has_alloced_mem = true;
  // assignment array holds the task assignment
  // reset to all -1

  hg_step3(costMatrix, sa, n);

  // now we've assigned as many values as we can
  // we'll check for a valid assignment
  int s4status;
  while (true){
    s4status = 0;
    memset(cc, 0, n);
    memset(cr, 0, n);
    do {
      s4status = hg_step4(sa, cc, cr, n, s4status);
    } while (s4status != 0);
    if (hg_statuscheck(sa, assignment, n)) break;
    hg_step5(costMatrix, cc, cr, n);
    hg_step3(costMatrix, sa, n);
    R_CheckUserInterrupt();
  }

  // if we made it here, we should have an optimal assignment
  free(sa);
  free(cc);
  free(cr);
  has_alloced_mem = false;
  return; 
}

void hg_step3(double *costMatrix, uint8_t *statusarray, int n){
  for(int i=0; i<n*n; i++) statusarray[i] = 0;

  bool foundzero, assigned;
  int idx;
  for(int i=0; i<n; i++){
    foundzero = false;
    for(int j=0; j<n; j++){
      idx = i*n+j;
      if (costMatrix[idx] == 0){
        // update bit flags
        statusarray[idx] = 3;
        if (!foundzero){
          // if haven't yet found a zero, try to assign
          assigned = false;
          for (int k=0; k<i; k++){
            if(statusarray[k*n+j] == 1){
              assigned = true;
              break;
            }
          }
          // if we this column hasn't been assigned, assign it
          if (!assigned){
            foundzero = true;
            statusarray[idx] = 1; 
          }
        }
      }
    }
  }
}

int hg_step4(uint8_t *statusarray, bool *coveredCol, bool *coveredRow, int n, int statuscode){
  // 0=none, 1=starred, 2=primed, 3=zero

  // first we cover all columns containing a starred zero
  // variables swapped so I always use j for columns and i for rows
  if (statuscode < 2){
    for (int j=0; j<n; j++){ // column iterator
      for (int i=0; i<n; i++){ // row iterator
        if(statusarray[i*n+j] == 1){
          coveredCol[j] = true;
          break;
        }
      }
    }
  }

  // Next, we find a non-covered zero and prime it
  for(int i=0; i<n; i++){
    if(coveredRow[i]) continue;
    for (int j=0; j<n; j++){
      if(coveredCol[j]) continue;
      if(statusarray[i*n+j] > 1){
        // we found an uncovered zero
        statusarray[i*n+j] = 2;

        // we then check if the zero is on the same row as a starred zero
        int starcol = 0;
        while(starcol<n && statusarray[i*n+starcol] != 1) starcol++;

        if(starcol < n){
          // we found a starred zero on the same row, 
          // so we uncover the row and cover the starred column
          coveredCol[starcol] = false;
          coveredRow[i] = true;
          // return 2 to skip covering columns on next iteration
          return 2;

        } else {
          // we didn't find a starred zero on the same row,
          // so we make a path alternating between primes and stars
          int starrow = i;
          int primecol = j;
          int tmp;
          while(true){
            // mark primed zero
            statusarray[starrow*n+primecol] += 2;

            // find a starred zero on same column (if we can't, break)
            tmp=0;
            while(tmp<n && statusarray[tmp*n+primecol] != 1) tmp++;
            if (tmp == n) break;
            statusarray[tmp*n+primecol] += 5;

            // find a primed zero on same row
            starrow=tmp;
            tmp=0;
            while(tmp<n && statusarray[starrow*n+tmp] != 2) tmp++;
            primecol=tmp;
          }

          // Once we've completed a path, we do the following:
          for (int r=0; r<n; r++){
            // uncovered all rows/columns
            coveredCol[r] = false;
            coveredRow[r] = false;
            for (int c=0; c<n; c++){
              // star primed zeros from path, unstar starred zeros
              // this will be either 6 (marked star) or 4 (marked prime)
              // leading to 3 (regular zero) or 1 (starred zero)
              if (statusarray[r*n+c] > 3) statusarray[r*n+c] -= 3;

              // unprime all primed zeros
              if (statusarray[r*n+c] == 2) statusarray[r*n+c]++;
            }
          }

          return 1;
        }
      }
    }
  }

  // if we made it here, we didn't find an uncovered zero
  return 0;
}

void hg_step5(double *costMatrix, bool *coveredCol, bool *coveredRow, int n){
  double minv = -1;

  for(int i=0; i<n; i++){
    if (coveredRow[i]) continue;
    for (int j=0; j<n; j++){
      if (coveredCol[j]) continue;
      minv = minv < 0 || costMatrix[i*n+j] < minv ? costMatrix[i*n+j] : minv;
    }
  }

  if (minv <= 0) return;
  // subtract this value from each uncovered row
  for(int i=0; i<n; i++){
    if (coveredRow[i]) continue;
    for (int j=0; j<n; j++){
      costMatrix[i*n+j] -= minv;
    }
  }

  // add this value to each covered column
  for(int j=0; j<n; j++){
    if (!coveredCol[j]) continue;
    for (int i=0; i<n; i++){
      costMatrix[i*n+j] += minv;
    }
  }

  return;
}

bool hg_statuscheck(uint8_t *statusarray, int *assignment, int n){
  for (int i=0; i<n; i++){
    assignment[i] = -1;
    for (int j=0; j<n; j++){
      if (statusarray[i*n+j] == 1)
        assignment[i] = j;
    }
    if (assignment[i] == -1) return false;
  }
  return true;
}



