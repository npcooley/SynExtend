#include "fortran_includes.h"
#include "RandomForest.h"
#include "SynExtend.h"

/*
 * Tree node initialization
 */
DTN *initNode(){
  DTN *node = malloc(sizeof(DTN));
  node->threshold = 0;
  node->index = -1;
  node->gini_gain = 0;
  node->left = NULL;
  node->right = NULL;
  return node;
}

/*
 * R-exposed functions
 */
SEXP R_learn_tree(SEXP DATA, SEXP NROWS, SEXP NCOLS, SEXP RESPONSE, SEXP NCLASSES,
  SEXP TO_CHECK, SEXP MAX_DEPTH, SEXP MIN_NODESIZE, SEXP ISCLASSIFICATION){
  int isclassif = LOGICAL(ISCLASSIFICATION)[0];

  // array input
  double *data = REAL(DATA);

  // variable inputs
  int nrows = INTEGER(NROWS)[0];
  int ncols = INTEGER(NCOLS)[0];
  int nclasses = INTEGER(NCLASSES)[0];
  int num_to_check = INTEGER(TO_CHECK)[0];
  int max_depth = INTEGER(MAX_DEPTH)[0];
  int min_nodesize = INTEGER(MIN_NODESIZE)[0];

  // internal vars
  DTN *head = initNode();

  // helper function will destroy data and class_response, so duplicate them first
  double *dup_data = malloc(sizeof(double)*nrows*ncols);

  // these do not need to be free'd -- will be free'd in the helper function
  //Rprintf("Duplicating memory.\n");
  dup_data = memcpy(dup_data, data, sizeof(double)*nrows*ncols);
  if(isclassif){
    int *response = INTEGER(RESPONSE);
    int *dup_response = malloc(sizeof(int)*nrows);
    dup_response = memcpy(dup_response, response, sizeof(int)*nrows);
    learntreeclassif_helper(head, dup_data, dup_response, nrows, ncols, nclasses,
                            num_to_check, 0, max_depth, min_nodesize);
  } else {
    double *response = REAL(RESPONSE);
    double *dup_response = malloc(sizeof(double)*nrows);
    dup_response = memcpy(dup_response, response, sizeof(double)*nrows);
    learntreeregress_helper(head, dup_data, dup_response, nrows, ncols,
                            num_to_check, 0, max_depth, min_nodesize);
  }

  // now we should have our entire tree created, and our duplicated arrays destroyed.

  // these objects will be allocated in `export_internal_tree`
  int *indices = NULL;
  double *thresholds = NULL, *gini_gain=NULL;
  int l = 0;
  export_internal_tree(head, &indices, &thresholds, &gini_gain, &l);

  // This is one option, I'm instead just going to register the external
  // pointer right away and return it, since I think that's easier.
  // Avoids a double call, and most people will predict right after
  // training anyway.

  // Read values back into R
  SEXP R_retval = PROTECT(allocVector(VECSXP, 4));
  SEXP R_indices = PROTECT(allocVector(INTSXP, l));
  SEXP R_thresholds = PROTECT(allocVector(REALSXP, l));
  SEXP R_gini = PROTECT(allocVector(REALSXP, l));

  memcpy(INTEGER(R_indices), indices, sizeof(int)*l);
  memcpy(REAL(R_thresholds), thresholds, sizeof(double)*l);
  memcpy(REAL(R_gini), gini_gain, sizeof(double)*l);
  free(indices);
  free(thresholds);

  SET_VECTOR_ELT(R_retval, 1, R_indices);
  SET_VECTOR_ELT(R_retval, 2, R_thresholds);
  SET_VECTOR_ELT(R_retval, 3, R_gini);
  UNPROTECT(3);

  // register the external pointer and then return
  SEXP R_ptr = PROTECT(R_MakeExternalPtr(head, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(R_ptr, (R_CFinalizer_t) R_TreeFinalizer, TRUE);
  SET_VECTOR_ELT(R_retval, 0, R_ptr);
  UNPROTECT(2);

  return(R_retval);
}

SEXP R_get_treeptr(SEXP VolatilePtr, SEXP INDICES, SEXP THRESHOLDS, SEXP GINIS){
  // We're going to lazily evaluate these trees
  // if they exist, just return the external pointer
  // note that it seems NULL can be treated as an external pointer address for whatever reason
  if(VolatilePtr != R_NilValue && R_ExternalPtrAddr(VolatilePtr)) return(VolatilePtr);
  int madePtr = 0;

  // otherwise, create the tree and return an external pointer to it
  DTN *tree = bfs_q2tree(INTEGER(INDICES), REAL(THRESHOLDS), REAL(GINIS), LENGTH(INDICES));
  if(VolatilePtr == R_NilValue){
    VolatilePtr = PROTECT(R_MakeExternalPtr(tree, R_NilValue, R_NilValue));
    madePtr = 1;
  } else {
    R_SetExternalPtrAddr(VolatilePtr, tree);
  }
  R_RegisterCFinalizerEx(VolatilePtr, (R_CFinalizer_t) R_TreeFinalizer, TRUE);
  if(madePtr) UNPROTECT(1);
  return VolatilePtr;
}

SEXP R_rfpredict(SEXP RF_Obj, SEXP DATA, SEXP L, SEXP NENTRIES){
  // assume we transpose the input before it's passed in so that each column is one entry
  // pointer needs to be reprotected because it may not be protected after R_get_treeptr
  SEXP R_ptr = PROTECT(R_get_treeptr(VECTOR_ELT(RF_Obj, 0), VECTOR_ELT(RF_Obj, 1), VECTOR_ELT(RF_Obj, 2), VECTOR_ELT(RF_Obj, 3)));
  SET_VECTOR_ELT(RF_Obj, 0, R_ptr);
  DTN *tree = (DTN *) R_ExternalPtrAddr(R_ptr);

  int nentries = INTEGER(NENTRIES)[0];
  int len = INTEGER(L)[0];
  double *data = REAL(DATA);
  SEXP R_return = PROTECT(allocVector(REALSXP, nentries));
  double *retval = REAL(R_return);

  for(int i=0; i<nentries; i++)
    retval[i] = predict_for_input(tree, &data[i*len]);

  UNPROTECT(2);
  return(R_return);
}


/*
 * Import/Export trees between C and R
 */
DTN *bfs_q2tree(int *indices, double *thresholds, double *gini, int length){
  queue *q = malloc(sizeof(queue));
  queue *end = q;
  queue *tmp_q = q;
  DTN *tmp, *head;

  head = initNode();
  q->ptr = head;
  q->next = NULL;
  int i=0, cur_ind;

  while(q && i<length){
    // load value into queue
    cur_ind = indices[i];
    tmp = q->ptr;
    tmp->threshold = thresholds[i];
    tmp->gini_gain = gini[i];
    tmp->index = cur_ind;
    if(cur_ind > -1){
      // add both children of the node into the queue
      end->next = malloc(sizeof(queue));
      end = end->next;
      tmp->left = initNode();
      end->ptr = tmp->left;
      end->next = malloc(sizeof(queue));
      end=end->next;
      tmp->right = initNode();
      end->ptr = tmp->right;
      end->next = NULL;
    }

    i++;
    q = q->next;
  }

  // free the entire queue
  while(tmp_q){
    q = tmp_q;
    tmp_q = tmp_q->next;
    free(q);
  }
  return head;
}

void export_internal_tree(DTN *tree, int **indices, double **thresholds, double **gini_gain, int *outlength){
  // Notable here: `indices` and `thresholds` are NOT ALLOCATED
  // They should be allocated in this function and then returned
  int outlen = 0;
  queue *q = malloc(sizeof(queue));
  queue *end = q;
  queue *tmp_q = q;

  q->ptr = tree;
  q->next = NULL;

  // create a BFS queue of the tree
  while(tmp_q){
    if(tmp_q->ptr->index != -1){
      end->next = malloc(sizeof(queue));
      end = end->next;
      end->ptr = tmp_q->ptr->left;
      end->next = malloc(sizeof(queue));
      end = end->next;
      end->ptr = tmp_q->ptr->right;
      end->next = NULL;
    }
    tmp_q = tmp_q->next;
    outlen++;
  }

  // allocate storage to be returned
  *indices = malloc(sizeof(int)*outlen);
  *thresholds = malloc(sizeof(double)*outlen);
  *gini_gain = malloc(sizeof(double)*outlen);

  // write values into storage (will be freed later)
  *outlength = outlen;
  for(int i=0; i<outlen; i++){
    (*indices)[i] = q->ptr->index;
    (*thresholds)[i] = q->ptr->threshold;
    (*gini_gain)[i] = q->ptr->gini_gain;
    tmp_q = q;
    q = q->next;
    free(tmp_q);
  }
  return;
}

/*
 * Other Helper Functions
 */
double predict_for_input(DTN *tree, double *data){
  DTN *tmp=tree;

  while(tmp->index != -1){
    if(data[tmp->index] <= tmp->threshold)
      tmp = tmp->left;
    else
      tmp = tmp->right;
  }

  return(tmp->threshold);
}

void learntreeclassif_helper(DTN *node, double *data, int *class_response,
                              int nrows, int ncols, int nclasses, int num_to_check,
                              int cur_depth, int max_depth, int min_nodesize){
  R_CheckUserInterrupt();
  // IMPORTANT: *data is assumed to malloc'd elsewhere
  // *data WILL BE FREE'd IN THIS FUNCTION
  // the same is true of *class_response
  // Calling function should duplicate any memory that cannot be destroyed

  // First we'll check to see if all the entries are the same class
  int curclass = nrows==0 ? -1 : class_response[0];
  int foundDifferent=0;
  for(int i=0; i<nrows; i++){
    if(curclass != class_response[i]){
      foundDifferent=1;
      break;
    }
  }
  if(!foundDifferent) cur_depth = max_depth;

  // if already at max_depth or we have fewer observations than minimum node size,
  // just assign the most prominent class and return
  if(cur_depth == max_depth || nrows <= min_nodesize){
    // these should be preinitialized but it doesn't hurt to be safe
    node->index=-1;
    node->threshold = 1; // default value is the first class
    node->left = NULL;
    node->right = NULL;
    node->gini_gain = 0;

    if(!foundDifferent){
      // find the most prominent class
      int *classcounts = calloc(nclasses, sizeof(int));
      for(int i=0; i<nrows; i++)
        classcounts[class_response[i]-1]++;

      for(int i=1; i<nclasses; i++)
        if(classcounts[i] > classcounts[(int)(node->threshold-1)])
          node->threshold = (double)(i+1);

      free(classcounts);
    } else {
      // if we already know they're all the same, no need to check all again
      node->threshold = (double)curclass;
    }

    free(class_response);
    free(data);
    return;
  }

  // otherwise we need to split into nodes
  split_decision_node_classif(node, data, class_response,
                              nrows, ncols, nclasses, num_to_check);

  double splitpoint = node->threshold;
  int ind = node->index;
  int nrow_left = 0, nrow_right=0;
  double *v = &data[nrows*ind];
  if(ind == -1){
    node->left = NULL;
    node->right = NULL;
    free(data);
    free(class_response);
    return;
  }

  // How big do we need the new arrays to be?
  for(int i=0; i<nrows; i++){
    if(v[i] <= splitpoint)
      nrow_left++;
    else
      nrow_right++;
  }

  double *left_data = malloc(sizeof(double) * nrow_left*ncols);
  double *right_data = malloc(sizeof(double) * nrow_right*ncols);
  int *left_class = malloc(sizeof(int) * nrow_left);
  int *right_class = malloc(sizeof(int) * nrow_right);
  int ctr_l=0, ctr_r=0;
  for(int i=0; i<nrows*ncols; i++){
    if(v[i%nrows] <= splitpoint){
      left_data[ctr_l] = data[i];
      if(ctr_l < nrow_left)
        left_class[ctr_l] = class_response[i%nrows];
      ctr_l++;
    } else {
      right_data[ctr_r] = data[i];
      if(ctr_r < nrow_right)
        right_class[ctr_r] = class_response[i%nrows];
      ctr_r++;
    }
  }

  // FREEING INPUT DATA/CLASS_RESPONSE
  free(data);
  free(class_response);

  DTN *left_node = initNode();
  DTN *right_node = initNode();

  // left node
  learntreeclassif_helper(left_node, left_data, left_class, nrow_left,
                          ncols, nclasses, num_to_check, cur_depth+1,
                          max_depth, min_nodesize);
  // right node
  learntreeclassif_helper(right_node, right_data, right_class, nrow_right,
                          ncols, nclasses, num_to_check, cur_depth+1,
                          max_depth, min_nodesize);

  node->left = left_node;
  node->right = right_node;

  return;
}

void split_decision_node_classif(DTN *node, double *data, int *class_response,
                                  int nrows, int ncols, int nclass, int num_to_check){
  // data should always be a numeric
  // response should be an int ranging from 1:n
  // nclass, num_to_check are constant throughout execution of the program

  // data will be a matrix stored by column (first nrows entries are col1, second are col2, etc.)
  // we'll just assume that all the preprocessing is done in R, no need to fiddle with that here
  // processing the SEXPs will be done separately so we can repeatedly call this internally

  // setting up a random sample of ints
  int *cols = malloc(sizeof(int) * ncols);
  for(int i=0; i<ncols; i++) cols[i] = i;
  int choice, tmp;

  // shuffle the columns
  GetRNGstate();
  for(int i=ncols-1; i>0; i--){
    choice = floor(unif_rand()*i);
    tmp = cols[choice];
    cols[choice] = cols[i];
    cols[i] = tmp;
  }
  PutRNGstate();

  double *results = malloc(sizeof(double) * num_to_check);
  double *gini_gain = malloc(sizeof(double) * num_to_check);
  double curmax = 0;
  choice = -1;
  for(int i=0; i<num_to_check; i++){
    F77_CALL(find_gini_split)(&data[nrows*cols[i]], class_response, &nrows, &nclass, &results[i], &gini_gain[i]);
    if(gini_gain[i] >= curmax){
      // using geq in case we find a split that also leaves gini unchanged
      // a split with equal gini is better than just terminating (probably?)
      choice = i;
      curmax = gini_gain[i];
    }
  }
  if(choice == -1){
    // assign the most prominent class here
    int *class_counts = calloc(nclass, sizeof(int));
    int whichmax = 0, countmax = -1;
    for(int i=0; i<nrows; i++){
      tmp = class_response[i]-1;
      class_counts[tmp]++;
      if(class_counts[tmp] > countmax){
        countmax = class_counts[tmp];
        whichmax = tmp+1;
      }
    }
    free(class_counts);
    node->index = -1;
    node->threshold = whichmax;
    node->gini_gain = 0.0;
  } else {
    node->threshold = results[choice];
    node->index = cols[choice];
    node->gini_gain = curmax;
  }

  free(results);
  free(gini_gain);
  free(cols);

  return;
}

void learntreeregress_helper(DTN *node, double *data, double *response,
                              int nrows, int ncols, int num_to_check,
                              int cur_depth, int max_depth, int min_nodesize){
  // this is a modification of learntreeclassif_helper
  R_CheckUserInterrupt();
  // IMPORTANT: *data is assumed to malloc'd elsewhere
  // *data WILL BE FREE'd IN THIS FUNCTION
  // the same is true of *response
  // Calling function should duplicate any memory that cannot be destroyed

  // First we'll check to see if all the entries are the same value
  double curval = nrows==0 ? -1 : response[0];
  int foundDifferent=0;
  for(int i=0; i<nrows; i++){
    if(curval != response[i]){
      foundDifferent=1;
      break;
    }
  }
  if(!foundDifferent) cur_depth = max_depth;

  // if already at max_depth or we have fewer observations than minimum node size,
  // just assign the most prominent class and return
  if(cur_depth == max_depth || nrows <= min_nodesize){
    // these should be preinitialized but it doesn't hurt to be safe
    node->index=-1;
    node->threshold = 1; // default value is the first class
    node->left = NULL;
    node->right = NULL;
    node->gini_gain = 0;

    if(!foundDifferent){
      // assign the average value
      curval = 0;
      for(int i=0; i<nrows; i++)
        curval += response[i] / nrows;
      node->threshold = curval;

    } else {
      // if we already know they're all the same, no need to check all again
      node->threshold = curval;
    }

    free(response);
    free(data);
    return;
  }

  // otherwise we need to split into nodes
  split_decision_node_regress(node, data, response,
                              nrows, ncols, num_to_check);

  double splitpoint = node->threshold;
  int ind = node->index;
  int nrow_left = 0, nrow_right=0;
  double *v = &data[nrows*ind];
  if(ind == -1){
    node->left = NULL;
    node->right = NULL;
    free(data);
    free(response);
    return;
  }

  // How big do we need the new arrays to be?
  for(int i=0; i<nrows; i++){
    if(v[i] <= splitpoint)
      nrow_left++;
    else
      nrow_right++;
  }

  double *left_data = malloc(sizeof(double) * nrow_left*ncols);
  double *right_data = malloc(sizeof(double) * nrow_right*ncols);
  double *left_response = malloc(sizeof(double) * nrow_left);
  double *right_response = malloc(sizeof(double) * nrow_right);
  int ctr_l=0, ctr_r=0;
  for(int i=0; i<nrows*ncols; i++){
    if(v[i%nrows] <= splitpoint){
      left_data[ctr_l] = data[i];
      if(ctr_l < nrow_left)
        left_response[ctr_l] = response[i%nrows];
      ctr_l++;
    } else {
      right_data[ctr_r] = data[i];
      if(ctr_r < nrow_right)
        right_response[ctr_r] = response[i%nrows];
      ctr_r++;
    }
  }

  // FREEING INPUT DATA/CLASS_RESPONSE
  free(data);
  free(response);

  DTN *left_node = initNode();
  DTN *right_node = initNode();

  // left node
  learntreeregress_helper(left_node, left_data, left_response, nrow_left,
                          ncols, num_to_check, cur_depth+1,
                          max_depth, min_nodesize);
  // right node
  learntreeregress_helper(right_node, right_data, right_response, nrow_right,
                          ncols, num_to_check, cur_depth+1,
                          max_depth, min_nodesize);

  node->left = left_node;
  node->right = right_node;

  return;
}

void split_decision_node_regress(DTN *node, double *data, double *response,
                                  int nrows, int ncols, int num_to_check){
  // data should always be a numeric
  // response should be a float
  // nclass, num_to_check are constant throughout execution of the program

  // data will be a matrix stored by column (first nrows entries are col1, second are col2, etc.)
  // we'll just assume that all the preprocessing is done in R, no need to fiddle with that here
  // processing the SEXPs will be done separately so we can repeatedly call this internally

  // split should be determined by minimizing sum of square error

  // setting up a random sample of ints
  int *cols = malloc(sizeof(int) * ncols);
  for(int i=0; i<ncols; i++) cols[i] = i;
  int choice, tmp;

  // shuffle the columns
  GetRNGstate();
  for(int i=ncols-1; i>0; i--){
    choice = floor(unif_rand()*i);
    tmp = cols[choice];
    cols[choice] = cols[i];
    cols[i] = tmp;
  }
  PutRNGstate();

  // here 'sse' stores the decrease in sse from before to after the split
  double *results = malloc(sizeof(double) * num_to_check);
  double *sse = malloc(sizeof(double) * num_to_check);
  double curmax = 0;
  choice = -1;
  for(int i=0; i<num_to_check; i++){
    F77_CALL(find_sse_split)(&data[nrows*cols[i]], response, &nrows, &results[i], &sse[i]);
    if(sse[i] >= curmax){
      // maximizing the decrease in sse
      // geq for same reasons as the classification case
      choice = i;
      curmax = sse[i];
    }
  }
  if(choice == -1){
    // assign mean value
    curmax = 0;
    for(int i=0; i<nrows; i++)
      curmax += response[i] / nrows;
    node->index = -1;
    node->threshold = curmax;
    node->gini_gain = 0.0;
  } else {
    node->threshold = results[choice];
    node->index = cols[choice];
    node->gini_gain = curmax;
  }

  free(results);
  free(sse);
  free(cols);

  return;
}

/*
 * Cleanup Functions
 */
void R_TreeFinalizer(SEXP TreePointer){
  if (!R_ExternalPtrAddr(TreePointer)) return;
  DTN *head = (DTN *) R_ExternalPtrAddr(TreePointer);
  freeDecisionTree(head);
  R_ClearExternalPtr(TreePointer);
  return;
}

void freeDecisionTree(DTN *tree){
  if(!tree) return;
  freeDecisionTree(tree->left);
  freeDecisionTree(tree->right);
  free(tree);
}