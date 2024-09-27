#ifndef SYNEXTEND_H
#define SYNEXTEND_H

// R_init_DECIPHER.c
void R_init_SynExtend(DllInfo *info);

/**** calcMIR2C.c ****/
void cleanupFxn(void);
SEXP calcMIcVec(SEXP V1, SEXP V2, SEXP UV, SEXP PSEUDOCOUNT);
SEXP calcMIVec(SEXP V1, SEXP V2, SEXP LEN);
SEXP trimCovar(SEXP fMAT, SEXP fSP, SEXP sSP, SEXP NV, SEXP NR);

/**** CDend.c ****/
// utility functions
SEXP initCDend(SEXP dend);
SEXP hashString(SEXP label);
SEXP printTree(SEXP tnPtr);
SEXP calcAllTreeLengths(SEXP tnPtr);

// scoring functions for EvoWeaver
SEXP calcGainLoss(SEXP tnPtr, SEXP occVec, SEXP convertToGL);
SEXP calcScoreGL(SEXP tnPtr, SEXP glv1, SEXP glv2);
SEXP calcScoreJaccard(SEXP ov1, SEXP ov2, SEXP NN);
SEXP calcScoreHamming(SEXP ov1, SEXP ov2, SEXP NN, SEXP norm);
SEXP cladeCollapsePA(SEXP tnPtr, SEXP ANCESTRAL_STATES);

// Tree distance
SEXP GRFInfo(SEXP tnPtr1, SEXP tnPtr2, SEXP allLabels, SEXP shouldUseJRF, SEXP JRFExp);
SEXP RFDist(SEXP tnPtr1, SEXP tnPtr2, SEXP allLabels);
SEXP KFDist(SEXP tnPtr1, SEXP tnPtr2, SEXP allLabels);

// D value calculation
SEXP calcDValue(SEXP tnPtr, SEXP occVec);
SEXP calcDRandValue(SEXP tnPtr, SEXP allLabels, SEXP numP, SEXP iterNum);
SEXP calcDBrownValue(SEXP tnPtr, SEXP allLabels, SEXP iterNum, SEXP SD, SEXP START, SEXP THRESH);

/**** XORRand.c ****/
SEXP pseudoRandomSample(SEXP N);
SEXP randomProjection(SEXP VEC, SEXP NONZERO, SEXP N, SEXP OUTDIM, SEXP NTHREADS);
SEXP seededPseudoRandomSample(SEXP N, SEXP SEED);
SEXP se_cophenetic(SEXP Index1, SEXP Index2, SEXP N, SEXP D, SEXP H);

/**** MoranI.c ****/
SEXP MoransI(SEXP VALS, SEXP DIST, SEXP DIM);

/**** CShuffle.c ****/
void shuffleRInt(int *v, int *l);
void shuffleRRepl(int *v, int *l);


/**** NucleotideCounts.c ****/
SEXP StringToNVDT(SEXP DNASTRING, SEXP REMOVEGAPS, SEXP EXTENDED, SEXP USEDNA);
SEXP fastPearsonC(SEXP V1, SEXP V2);
SEXP MIForSequenceSets(SEXP M1, SEXP M2, SEXP NSEQS, SEXP U1, SEXP U2, SEXP BASE, SEXP NTHREADS);

/**** HungarianAlgo.c ****/
SEXP HungarianAssignment(SEXP MATVEC, SEXP DIM);
void hungarianCleanup(void);

/**** Dendrapply ****/
void free_dendrapply_list(void);
SEXP do_dendrapply(SEXP tree, SEXP fn, SEXP env, SEXP order);

/**** OnDiskLP.c ****/
SEXP R_LPOOM_cluster(SEXP FILENAME, SEXP NUM_EFILES, SEXP TABNAME,
                    SEXP TEMPTABNAME, SEXP QFILE, SEXP OUTDIR, SEXP OUTFILE,
                    SEXP SEPS, SEXP CTR, SEXP ITER, SEXP VERBOSE,
                    SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS, SEXP IGNORE_WEIGHTS, SEXP NORMALIZE_WEIGHTS,
                    SEXP CONSENSUS_WEIGHTS, SEXP INFLATION_POW, SEXP SHUFFLE_QUEUES);

/**** RandomForest.c ****/
SEXP R_learn_tree(SEXP DATA, SEXP NROWS, SEXP NCOLS, SEXP RESPONSE,
                          SEXP NCLASSES, SEXP TO_CHECK, SEXP MAX_DEPTH,
                          SEXP MIN_NODESIZE, SEXP ISCLASSIFICATION);
SEXP R_get_treeptr(SEXP VolatilePtr, SEXP INDICES, SEXP THRESHOLDS, SEXP GINIS);
SEXP R_rfpredict(SEXP RF_Obj, SEXP DATA, SEXP L, SEXP NENTRIES);

#endif
