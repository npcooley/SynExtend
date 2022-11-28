// R_init_DECIPHER.c
void R_init_SynExtend(DllInfo *info);

/**** calcMIR2C.c ****/
void cleanupFxn();
SEXP calcMIcVec(SEXP V1, SEXP V2, SEXP UV, SEXP PSEUDOCOUNT);
SEXP trimCovar(SEXP fMAT, SEXP fSP, SEXP sSP, SEXP NV, SEXP NR);

/**** CDend.c ****/
// utility functions
SEXP initCDend(SEXP dend);
SEXP hashString(SEXP label);
SEXP printTree(SEXP tnPtr);

// scoring functions for ProtWeaver
SEXP calcGainLoss(SEXP tnPtr, SEXP occVec, SEXP convertToGL);
SEXP calcScoreGL(SEXP tnPtr, SEXP glv1, SEXP glv2);
SEXP calcScoreJaccard(SEXP ov1, SEXP ov2, SEXP NN);
SEXP calcScoreHamming(SEXP ov1, SEXP ov2, SEXP NN, SEXP norm);

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
SEXP randomProjection(SEXP VEC, SEXP NONZERO, SEXP N, SEXP OUTDIM);
SEXP seededPseudoRandomSample(SEXP N, SEXP SEED);