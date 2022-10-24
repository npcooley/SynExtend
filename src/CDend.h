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
void findNodeScores(treeNode* curNode, int* v1, int* v2, double* scores, treeNode* head, bool isHead);
treeNode* findNextNode(treeNode *curNode, int *v, int *selfv, bool isCur);

/* Fitch Parsimony */
void resetTree(treeNode* node, int val);
void fitchUp(treeNode* node, unsigned int* hashMap, int hashMapLen, int* PAvec);
void fitchDown(treeNode* node, int parentVal, int* PAvec);
void fitchRecon(int* PAvec, int len, int defaultVal);
void convertGL(treeNode* node, bool curVal, int* PAvec);
int populateVector(treeNode* node, int* container, int idx);

/* internal functions */
static inline int getNumNodes(treeNode* node, int n);
static inline int labelTreePostorder(treeNode* node, int n);
static inline void checkPtrExists(SEXP tnPtr);
static void FreeTree(SEXP tnPtr);
static void printHelper(treeNode* node, int depth);
static void CleanupTree(treeNode* head);

SEXP TREEHT, TREEMEM, TREELAB, TREELF;
