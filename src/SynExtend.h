// R_init_DECIPHER.c

void R_init_SynExtend(DllInfo *info);

// calcMIR2C.c
void cleanupFxn();
SEXP calcMIcVec(SEXP V1, SEXP V2, SEXP UV, SEXP PSEUDOCOUNT);
SEXP trimCovar(SEXP fMAT, SEXP fSP, SEXP sSP, SEXP NV, SEXP NR);


// CDend.c
SEXP initCDend(SEXP dend);
SEXP hashString(SEXP label);
SEXP printTree(SEXP tnPtr);
SEXP calcGainLoss(SEXP tnPtr, SEXP occVec, SEXP convertToGL);
SEXP calcScoreGL(SEXP tnPtr, SEXP glv1, SEXP glv2);
SEXP calcScoreJaccard(SEXP ov1, SEXP ov2, SEXP NN);
SEXP calcScoreHamming(SEXP ov1, SEXP ov2, SEXP NN, SEXP norm);