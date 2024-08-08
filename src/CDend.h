#ifndef CDEND_H
#define CDEND_H

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

typedef unsigned long ulong;

typedef struct tree {
  double height;
  double inc_len;
  int value;
  int members;
  unsigned int label;
  struct tree *left, *right;
} treeNode;

treeNode *allocTreeNode(double h, int v, int m, unsigned int l, double parentheight);
treeNode *convertRDend(SEXP dend, double parentheight);

/* helpers for external functions */
const char *convertRChar(SEXP chars);
unsigned int hashLabel(SEXP label);
void findNodeScores(treeNode* curNode, int* v1, int* v2, double* scores, treeNode* head, bool isHead);
treeNode* findNextNode(treeNode *curNode, int *v, int *selfv, bool isCur);
void cladeCollapseHelper(treeNode* node, int *v);
void treeLengthsHelper(treeNode* head, double *v);

/* D value on Trees */
void calcSisterClades(treeNode *node, unsigned int *pmap, int pmaplen, double *scoreArr);
double scoreSisterClades(treeNode *node, double *scores);
void propBrownianEvo(treeNode *node, double *scores, double curval, double sd);
void findMapping(treeNode *node, int *mapping, unsigned int *hashvals, int lenHash);

/* Tree Distance */
void internalPartitionMap(treeNode *node, bool **pSets, unsigned int *hvs, int lh, int rootv);
int reallocPartitionMap(bool **pSets, int lh, int plen);
double scorePMs(bool **pm1, bool **pm2, int pm1l, int pm2l, int lh);
double calcEntropy(bool **pm, int lh, int pml);
ulong RFHashMap(treeNode *node, ulong *htable, ulong *keys, unsigned int *hvs, int lh, int rootv);
ulong KFHashMap(treeNode *node, ulong *htable, double *dists, ulong *keys, unsigned int *hvs, int lh, int rootv);
double scoreJaccardRFDist(bool **pm1, bool **pm2, int pm1l, int pm2l, int lh, double expv);
double calcJaccardPairingScore(bool *v1, bool *v2, int lh, double expv);

/* Fitch Parsimony */
void resetTree(treeNode* node, int val);
void fitchUp(treeNode* node, unsigned int* hashMap, int hashMapLen, int* PAvec);
void fitchDown(treeNode* node, int parentVal, int* PAvec);
void fitchRecon(int* PAvec, int len, int defaultVal);
void convertGL(treeNode* node, bool curVal, int* PAvec);
int populateVector(treeNode* node, int* container, int idx);

/* internal functions */
static inline double PclDist(double f, double s, double o, int tl){
  return o == 0 ? 0 : ((o/tl) * (log2((tl * o) / (s * f))));
}
static inline int getNumNodes(treeNode* node, int n);
static inline int labelTreePostorder(treeNode* node, int n);
static inline treeNode* checkPtrExists(SEXP tnPtr);
static void FreeTree(SEXP tnPtr);
static void printHelper(treeNode* node, int depth);
static void CleanupTree(treeNode* head);

SEXP TREEHT, TREEMEM, TREELAB, TREELF;

#endif