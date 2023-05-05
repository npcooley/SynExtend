/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

// for math functions
#include <math.h>

// SynExtend header file
#include "SynExtend.h"

// Other header files for .C Routines
#include "SEutils.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define C_DEF(name, n)  {#name, (DL_FUNC) &name, n}

/*
 * -- REGISTRATION OF THE .Call ENTRY POINTS ---
 */
static const R_CallMethodDef callMethods[] = { // method call, pointer, num args
  CALLDEF(calcMIcVec, 4),
  CALLDEF(calcMIVec, 3),
  CALLDEF(trimCovar, 5),
  CALLDEF(initCDend, 1),
  CALLDEF(hashString, 1),
  CALLDEF(calcGainLoss, 3),
  CALLDEF(calcScoreGL, 3), 
  CALLDEF(calcScoreJaccard, 3),
  CALLDEF(calcScoreHamming, 4),
  CALLDEF(printTree, 1),
  CALLDEF(GRFInfo, 5),
  CALLDEF(RFDist, 3),
  CALLDEF(KFDist, 3),
  CALLDEF(calcDValue, 2),
  CALLDEF(calcDRandValue, 4),
  CALLDEF(calcDBrownValue, 6),
  CALLDEF(pseudoRandomSample, 1),
  CALLDEF(randomProjection, 5),
  CALLDEF(seededPseudoRandomSample, 2),
  CALLDEF(MoransI, 3),
  CALLDEF(StringToNVDT, 4),
  CALLDEF(rpdendrapply, 3),
  CALLDEF(HungarianAssignment, 2),
  CALLDEF(fastPearsonC, 2),
  CALLDEF(do_dendrapply, 4),
  CALLDEF(se_cophenetic, 5),
  //CALLDEF(R_initNNptr, 6),
  //CALLDEF(R_PredictForInput, 2),
  //CALLDEF(R_UpdateWeights, 2),
  {NULL, NULL, 0}
};

/*
 * -- REGISTRATION OF THE .C ENTRY POINTS ---
 */
static const R_CMethodDef cMethods[] = {
  C_DEF(cleanupFxn, 0),
  C_DEF(shuffleRInt, 2),
  C_DEF(shuffleRRepl, 2),
  C_DEF(hungarianCleanup, 0),
  C_DEF(genCostMatrix, 7),
  C_DEF(free_dendrapply_list, 0),
  C_DEF(R_combineDistObj, 6),
  {NULL, NULL, 0}
};

void R_init_SynExtend(DllInfo *info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
