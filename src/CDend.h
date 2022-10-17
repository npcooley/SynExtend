#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "SynExtend.h"

typedef struct tree {
  double height;
  int value;
  int members;
  unsigned int label;
  struct tree *left, *right;
} treeNode;

treeNode *allocTreeNode(double h, int v, int m, unsigned int l);
treeNode *convertRDend(SEXP dend);

/* helpers for external functions */
const char *convertRChar(SEXP chars);
unsigned int hashLabel(SEXP label);
int findNodeScores(treeNode *curNode, int *v1, int *v2, double *scores, int ctr);
int findNextNode(treeNode *curNode, int val, int *v, int ctr);

/* Fitch Parsimony */
void resetTree(treeNode* node, int val);
void fitchUp(treeNode* node, unsigned int* hashMap, int hashMapLen);
void fitchDown(treeNode* node, int parentVal);
void fitchRecon(treeNode* node, int defaultVal);
void convertGL(treeNode* node, bool curVal);
int populateVector(treeNode* node, int *container, int idx);

/* internal functions */
static inline int getNumNodes(treeNode* node);
static inline void checkPtrExists(SEXP tnPtr);
static void FreeTree(SEXP tnPtr);
static void printHelper(treeNode* node, int depth);
static void CleanupTree(treeNode* head);

SEXP TREEHT, TREEMEM, TREELAB, TREELF;
treeNode *globalTreeNode;
bool globalIsSame;
