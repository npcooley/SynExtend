#ifndef RANDOMFOREST_INC
#define RANDOMFOREST_INC

#include <Rdefines.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// these structs will only be used in this file
// they should not be expected to persist

struct DTreeNode {
  // TODO: maybe save the node size so that we can prune better later
  struct DTreeNode *left;
  struct DTreeNode *right;
  double threshold;
  double gini_gain;
  int index;
};
typedef struct DTreeNode DTN;

struct DTNqueue{
  struct DTNqueue *next;
  DTN *ptr;
};
typedef struct DTNqueue queue;

DTN *initNode();

// internal functions
DTN *bfs_q2tree(int *indices, double *thresholds, double *gini, int length);

void R_TreeFinalizer(SEXP TreePointer);
void freeDecisionTree(DTN *tree);

void learntreeclassif_helper(DTN *node, double *data, int *class_response,
                              int nrows, int ncols, int nclass, int num_to_check,
                              int cur_depth, int max_depth, int min_nodesize);
void split_decision_node_classif(DTN *node, double *data, int *class_response,
                                  int nrows, int ncols, int nclass, int num_to_check);

void learntreeregress_helper(DTN *node, double *data, double *response,
                              int nrows, int ncols, int num_to_check,
                              int cur_depth, int max_depth, int min_nodesize);
void split_decision_node_regress(DTN *node, double *data, double *response,
                                  int nrows, int ncols, int num_to_check);

void export_internal_tree(DTN *tree, int **indices, double **thresholds, double **gini_gain, int *outlength);

double predict_for_input(DTN *tree, double *data);

/* testing stuff */
struct printQ {
  struct printQ *next;
  DTN *node;
  int index, depth;
  double threshold;
};
void printDecisionTree(DTN *tree);

#endif