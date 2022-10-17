#include "CDend.h"

/****** Callable Functions ******/
SEXP initCDend(SEXP dend){
  TREEHT = install("height");
  TREEMEM = install("members");
  TREELAB = install("label");
  TREELF = install("leaf");
  treeNode *head = convertRDend(dend);
  SEXP retval = PROTECT(R_MakeExternalPtr(head, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(retval, (R_CFinalizer_t) FreeTree, TRUE);
  UNPROTECT(1);
  return retval;
}

/* Prints the tree
Unneeded, but keeping in for future reference
(since I know I'll forget how to do this later)
SEXP printTree(SEXP tnPtr){
  checkPtrExists(tnPtr);
  treeNode *head = (treeNode *) R_ExternalPtrAddr(tnPtr);
  printHelper(head, 0);
  return R_NilValue;
}
*/

SEXP getTreeNodesCount(SEXP tnPtr){
  checkPtrExists(tnPtr);
  treeNode *head = (treeNode *) R_ExternalPtrAddr(tnPtr);
  SEXP retval = PROTECT(allocVector(INTSXP, 1));
  INTEGER(retval)[0] = getNumNodes(head);
  UNPROTECT(1);
  return retval;
}

SEXP hashString(SEXP label){
  SEXP ret = PROTECT(allocVector(INTSXP, 1));
  INTEGER(ret)[0] = hashLabel(STRING_ELT(label, 1));
  UNPROTECT(1);
  return ret;
}

SEXP calcGainLoss(SEXP tnPtr, SEXP occVec, SEXP convertToGL){
  checkPtrExists(tnPtr);
  treeNode *head = (treeNode *) R_ExternalPtrAddr(tnPtr);
  bool shouldConvert = LOGICAL(convertToGL)[0];
  int ovLen = LENGTH(occVec);
  unsigned int *presMap = malloc(sizeof(unsigned int) * ovLen);
  for (int i=0; i<ovLen; i++){
    presMap[i] = hashLabel(STRING_ELT(occVec, i));
  }

  resetTree(head, 2);
  fitchUp(head, presMap, ovLen);
  free(presMap);

  fitchDown(head, 0);
  fitchRecon(head, 0);
  if (shouldConvert)
    convertGL(head, head->value==1);

  int numNodes = getNumNodes(head);
  SEXP retvec = PROTECT(allocVector(INTSXP, numNodes));
  int *rvPtr = INTEGER(retvec);
  populateVector(head, rvPtr, 0);

  // commenting this out to save runtime
  // ALWAYS ASSUME TREE VALS HAVE NOT BEEN RESET BEFORE RUNNING OPERATIONS
  // resetTree(head, 2);
  UNPROTECT(1);
  return retvec;
}

SEXP calcScoreGL(SEXP tnPtr, SEXP glv1, SEXP glv2, SEXP NN){
  checkPtrExists(tnPtr);
  if (LENGTH(glv1) != LENGTH(glv2)) error("Gain/Loss vectors are different lengths!");
  treeNode *head = (treeNode *) R_ExternalPtrAddr(tnPtr);
  int *v1 = INTEGER(glv1);
  int *v2 = INTEGER(glv2);
  int numNodes = INTEGER(NN)[0];
  double *nodeScores = malloc(sizeof(double) * numNodes);
  findNodeScores(head, v1, v2, nodeScores, 0);

  
  double r = 0.0;
  for ( int i=0; i < numNodes; i++){
    if (nodeScores[i] != 0){
      r += 1 / sqrt(fabs(nodeScores[i]));
    }
  }

  free(nodeScores);
  globalTreeNode = NULL;
  SEXP retval = PROTECT(allocVector(REALSXP, 1));
  REAL(retval)[0] = r; 
  UNPROTECT(1);
  return retval;
}

SEXP calcScoreJaccard(SEXP ov1, SEXP ov2, SEXP NN){
  int v1len = INTEGER(NN)[0];
  int *v1 = INTEGER(ov1);
  int *v2 = INTEGER(ov2);

  double outval = 0.0;
  for (int i=0; i<v1len; i++)
    outval += v1[i] == v2[i];
  double denom = (double)(2*v1len) - outval;
  outval = outval / denom;

  SEXP retval = PROTECT(allocVector(REALSXP, 1));
  REAL(retval)[0] = outval;
  UNPROTECT(1);
  return retval;
}

/****** Gain / Loss Functions ******/
int findNodeScores(treeNode *curNode, int *v1, int *v2, double *scores, int ctr){
  scores[ctr] = 0.0;
  if (v1[ctr] != 0){
    globalTreeNode = NULL;
    findNextNode(curNode, v1[ctr], v2, ctr);
    if (globalTreeNode)
      scores[ctr] = (!globalIsSame * -2 + 1) * (fabs(curNode->height - globalTreeNode->height) + 1);
  }
  ctr++;
  if (curNode->left) ctr = findNodeScores(curNode->left, v1, v2, scores, ctr);
  if (curNode->right) ctr = findNodeScores(curNode->right, v1, v2, scores, ctr);

  return ctr;
}

int findNextNode(treeNode *curNode, int val, int *v, int ctr){
  if (v[ctr] != 0){
    if (!globalTreeNode || curNode->height < globalTreeNode->height){
      globalTreeNode = curNode;
      globalIsSame = v[ctr] == val;
    }
  }  
  ctr++;

  if (curNode->left) ctr = findNextNode(curNode->left, val, v, ctr);
  if (curNode->right) ctr = findNextNode(curNode->right, val, v, ctr);

  return ctr;
}


/****** Fitch Parsimony Functions ******/

void resetTree(treeNode* node, int val){
  node->value = val;
  if (node->left) resetTree(node->left, val);
  if (node->right) resetTree(node->right, val);
}

// Up Phase, starts with leaves and propagates up
void fitchUp(treeNode* node, unsigned int* hashMap, int hashMapLen){
  if (node->label != 0){
    // leaf node
    bool found = false;
    for (int i = 0; i<hashMapLen; i++){
      if (node->label == hashMap[i]){
        found = true;
        break;
      }
    }

    node->value = found ? 1 : 0;
    return;
  }

  // non-leaf node
  int lv=2, rv=2;
  if (node->left){
    fitchUp(node->left, hashMap, hashMapLen);
    lv = node->left->value;
  }

  if (node->right){
    fitchUp(node->right, hashMap, hashMapLen);
    rv = node->right->value;
  }

  // assign value based on children
  if (rv == 2 || lv == 2){
    node->value = rv==2 ? lv : rv;
  } else {
    node->value = rv==lv ? rv : 2;
  }

  return;
}

// Down phase, takes consensus between parent and children
void fitchDown(treeNode* node, int parentVal){
  if (node->label != 0) return;
  if (node->value == 2){
    int counts[3] = {0, 0, 0};
    counts[parentVal]++;

    int lv = node->left ? node->left->value : 2;
    int rv = node->right ? node->right->value : 2;
    counts[lv]++;
    counts[rv]++;

    int endval = 2;
    if (counts[2] != 3 && counts[0] != counts[1])
      endval = counts[1] > counts[0] ? 1 : 0;

    node->value = endval;
  }

  if (node->left) fitchDown(node->left, node->value);
  if (node->right) fitchDown(node->right, node->value);
}

// Reconciliation, any nodes we couldn't set we just set to defaultVal
void fitchRecon(treeNode* node, int defaultVal){
  if (node->value == 2) node->value = defaultVal;
  if (node->left) fitchRecon(node->left, defaultVal);
  if (node->right) fitchRecon(node->right, defaultVal);
  return;
}

void convertGL(treeNode* node, bool curVal){
  if ((node->value==1) ^ curVal){
    // if curVal was true and we've changed, it's a loss (else gain)
    node->value = curVal ? -1 : 1;
    curVal = !curVal;
  } else {
    node->value = 0;
  }

  if (node->left) convertGL(node->left, curVal);
  if (node->right) convertGL(node->right, curVal);
  return;
}

int populateVector(treeNode* node, int *container, int idx){
  int curIdx = idx;
  container[curIdx] = node->value;
  curIdx++;
  if (node->left) curIdx = populateVector(node->left, container, curIdx);
  if (node->right) curIdx = populateVector(node->right, container, curIdx);
  
  return curIdx;
}

/****** Allocation Functions ******/

treeNode *allocTreeNode(double h, int v, int m, unsigned int l){
  treeNode *newNode = Calloc(1, treeNode);
  newNode->left = NULL;
  newNode->right = NULL;
  newNode->height = h;
  newNode->value = v;
  newNode->members = m;
  newNode->label = l;
  return newNode;
}

/****** Conversion Functions ******/

// convert dendrogram to treeNode
treeNode *convertRDend(SEXP dend){
  double h = 0.0;
  int v = -1;
  int m = 1;
  unsigned int label = 0;

  if (!isNull(getAttrib(dend, TREEHT)))
    h = REAL(getAttrib(dend, TREEHT))[0];

  if (!isNull(getAttrib(dend, TREEMEM)))
    m = INTEGER(getAttrib(dend, TREEMEM))[0];

  if (!isNull(getAttrib(dend, TREELAB)))
      label = hashLabel(STRING_ELT(getAttrib(dend, TREELAB), 0));

  if (!isNull(getAttrib(dend, TREELF))){
    return allocTreeNode(h, v, m, label);
  }

  treeNode *cur = allocTreeNode(h, v, m, label);
  cur->left = convertRDend(VECTOR_ELT(dend, 0));
  cur->right = convertRDend(VECTOR_ELT(dend, 1));

  return cur;
}

unsigned int hashLabel(SEXP label){
  const char *s = Rf_translateCharUTF8(label);

  unsigned int result = 0x55555555;

  while (*s) { 
      result ^= *s++;
      result = ((result << 5) | (result >> (27)));
  }

  return result;
}

/****** Miscellaneous ******/

static inline int getNumNodes(treeNode* node){
  int numNodes = 1;
  if (node->left) numNodes += getNumNodes(node->left);
  if (node->right) numNodes += getNumNodes(node->right);
  return numNodes;
}

/* See comment on printTree
// helper function to print tree (in order)
static void printHelper(treeNode* node, int depth){
  if (!node) return;
  int vl = node->left ? node->left->value : 0;
  int vr = node->right ? node->right->value : 0;
  if (node->label != 0) Rprintf("LEAF %u (value %d)", node->label, node->value);
  else {
    Rprintf("%f (v: %d, l: %d, r: %d)", 
              node->height, node->value, vl, vr);
  }
  Rprintf("\n");
  if (node->left)
    printHelper(node->left, depth+1);
  if (node->right)
    printHelper(node->right, depth+1);

  return;
}
*/

static void FreeTree(SEXP tnPtr){
  if (!R_ExternalPtrAddr(tnPtr)) return;
  treeNode *head = (treeNode *) R_ExternalPtrAddr(tnPtr);
  CleanupTree(head);
  R_ClearExternalPtr(tnPtr);
  return;
}

static inline void checkPtrExists(SEXP tnPtr){
  if (!R_ExternalPtrAddr(tnPtr))
    error("External pointer no longer exists!");
  return;
}

// helper function to free memory
static void CleanupTree(treeNode* head){
  if (!head) return;
  CleanupTree(head->left);
  CleanupTree(head->right);
  Free(head);

  return;
}
