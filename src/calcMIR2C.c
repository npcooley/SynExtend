#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "SynExtend.h"

typedef struct node {
  double data;
  int i1;
  int i2;
  struct node *next;
} node;

static node *corrs = NULL;

void cleanupFxn(){
  node *ptr = corrs;
  // Free allocated memory in linked list
  while (corrs != NULL){
    ptr = corrs;
    corrs = corrs->next;
    free(ptr);
  }
}

node *makeNewNode(double data, int i1, int i2){
  node *newNode; 
  // very important to use calloc here we can free it later
  newNode = (node *)calloc(1, sizeof(node));
  if (newNode == NULL){
    error("Could not allocate storage for linked list.\n");
  }
  newNode->data = data;
  newNode->i1 = i1;
  newNode->i2 = i2;
  newNode->next = NULL;
  return(newNode);
}

// insert link in sorted order
// returns a bool corresponding to if we inserted
bool insertSorted(node **head, node *toInsert, int maxSize) {
  int ctr = 0;
  if ((*head) == NULL || (*head)->data >= toInsert->data){
    toInsert->next = *head;
    *head = toInsert;
    return(true);
  } else {
    node *temp = *head;
    while (temp->next != NULL && temp->next->data < toInsert->data){
      temp = temp->next;
      if (ctr == maxSize){
        return(false);
      }
      ctr += 1;
    }
    toInsert->next = temp->next;
    temp->next = toInsert;
    return(true);
  }
}

SEXP calcMIcVec(SEXP V1, SEXP V2, SEXP UV, SEXP PSEUDOCOUNT)
{
  SEXP ans = PROTECT(allocVector(REALSXP, 1));
  double *rans = REAL(ans);
  int l = length(V1);
  int uv = asInteger(UV);
  double pcount = REAL(PSEUDOCOUNT)[0];
  int *v1 = INTEGER(V1);
  int *v2 = INTEGER(V2);
  double *marg1 = (double*) S_alloc(uv, sizeof(double));
  double *marg2 = (double*) S_alloc(uv, sizeof(double));
  double *joint = (double*) S_alloc(uv*uv, sizeof(double));

  int val1, val2;
  for ( int i=0; i<l; i++ ){
    val1 = v1[i];
    val2 = v2[i];
    joint[val1*uv + val2] += 1;
    marg1[val1] += 1;
    marg2[val2] += 1;
  }

  double p1, p2, jp, outval;
  outval = 0;
  for ( int i=0; i<uv; i++ ){
    for ( int j=0; j<uv; j++ ){
      p1 = (marg1[i] + pcount) / (l + uv*pcount);
      p2 = (marg2[j] + pcount) / (l + uv*pcount);
      jp = (joint[i*uv + j] + pcount) / (l + pcount*uv*uv);
      if (p1 != 0 && p2 != 0 && jp != 0){
        outval += jp * log2(jp / (p1 * p2));
      }
    }
  }

  *rans = outval;
  UNPROTECT(1);
  return ans;
}

SEXP calcMIVec(SEXP V1, SEXP V2, SEXP LEN){//, SEXP PSEUDOCOUNT){
  // V1, V2 integer vectors
  // LEN is the length of V1, V2 (int)
  // UV is number of unique values (int)
  // PSEUDOCOUNT is pseudocounts to add (double)
  int *v1 = INTEGER(V1);
  int *v2 = INTEGER(V2);
  int l = LENGTH(V1);
  int uv = INTEGER(LEN)[0];
  //double pcount = REAL(PSEUDOCOUNT)[0];

  double mi = 0;
  double jointentropy = 0;

  double *pX, *pY, *pJoint;
  double div = ((double) 1) / ((double) l);
  pX = calloc(uv, sizeof(double));
  pY = calloc(uv, sizeof(double)); 
  pJoint = calloc(uv*uv, sizeof(double));
  int c1, c2;

  // calc individual probability dists
  for(int i=0; i<l; i++){
    c1 = v1[i]-1;
    c2 = v2[i]-1;
    pX[c1] += div;
    pY[c2] += div;
    pJoint[c1*uv+c2] += div;
  }

  // calc MI
  int idx;
  for(int i=0; i<uv; i++){
    for(int j=0; j<uv; j++){
      idx = i*uv+j;
      if(pJoint[idx] != 0){
        mi += pJoint[idx] * log2(pJoint[idx] / (pX[i]*pY[j]));
        jointentropy += pJoint[idx] * log2(pJoint[idx]);
      }
    }
  }

  mi /= (-1*jointentropy);
  SEXP out = PROTECT(allocVector(REALSXP, 1));
  REAL(out)[0] = mi;
  free(pX);
  free(pY);
  free(pJoint);
  UNPROTECT(1);
  return(out);
}


SEXP trimCovar(SEXP fMAT, SEXP fSP, SEXP sSP, SEXP NV, SEXP NR){
  int nv = asInteger(NV);
  int nr = asInteger(NR);
  int sp1l = length(fSP);
  int sp2l = length(sSP);

  int *firstSeqPos = INTEGER(fSP);
  int *secondSeqPos = INTEGER(sSP);
  double *fm = REAL(fMAT);
  int colv1, colv2;
  
  // Using a linked list for efficient insert
  corrs = NULL;
  int cv1, cv2;
  double p1, p2, score=0;
  bool success;
  for ( int i=0; i<sp1l; i++ ){
    cv1 = firstSeqPos[i];
    colv1 = (cv1 - 1) * nr;
    for ( int j=0; j<sp2l; j++ ){
      cv2 = secondSeqPos[j];
      colv2 = (cv2 - 1) * nr;

      score = 0;
      for ( int k=0; k<nr; k++){
        p1 = fm[colv1 + k];
        p2 = fm[colv2 + k];
        if (p1 != 0 && p2 != 0){
          score += p1 * log(p1 / p2);
        }
      }
      node *newNode = makeNewNode(score, cv1, cv2);
      success = insertSorted(&corrs, newNode, nv);
      // If we don't insert, free the associated memory
      if (!success && newNode != NULL){
        free(newNode);
        newNode = NULL;
      }
    }
    R_CheckUserInterrupt();
  }


  SEXP ans = PROTECT(allocVector(INTSXP, 2*nv));
  int *rans = INTEGER(ans); 
  node *ptr=corrs;

  for ( int i=0; i<nv; i++){
    rans[2*i] = ptr->i1;
    rans[2*i+1] = ptr->i2;
    ptr = ptr->next;
  }

  // Memory is freed after function terminates via R

  UNPROTECT(1);
  return(ans);
}








