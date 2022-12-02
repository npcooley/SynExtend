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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
// #include "Biostrings_interface.h"

// SynExtend header file
#include "SynExtend.h"

// Other header files for .C Routines
#include "SEutils.h"

/*
 * -- REGISTRATION OF THE .Call ENTRY POINTS ---
 */
static const R_CallMethodDef callMethods[] = { // method call, pointer, num args
  {"calcMIcVec", (DL_FUNC) &calcMIcVec, 4},
  {"trimCovar", (DL_FUNC) &trimCovar, 5},
  {"initCDend", (DL_FUNC) &initCDend, 1},
  {"hashString", (DL_FUNC) &hashString, 1},
  {"calcGainLoss", (DL_FUNC) &calcGainLoss, 3},
  {"calcScoreGL", (DL_FUNC) &calcScoreGL, 3}, 
  {"calcScoreJaccard", (DL_FUNC) &calcScoreJaccard, 3},
  {"calcScoreHamming", (DL_FUNC) &calcScoreHamming, 4},
  {"printTree", (DL_FUNC) &printTree, 1},
  {"GRFInfo", (DL_FUNC) &GRFInfo, 5},
  {"RFDist", (DL_FUNC) &RFDist, 3},
  {"KFDist", (DL_FUNC) &KFDist, 3},
  {"calcDValue", (DL_FUNC) &calcDValue, 2},
  {"calcDRandValue", (DL_FUNC) &calcDRandValue, 4},
  {"calcDBrownValue", (DL_FUNC) &calcDBrownValue, 6},
  {"pseudoRandomSample", (DL_FUNC) &pseudoRandomSample, 1},
  {"randomProjection", (DL_FUNC) &randomProjection, 4},
  {"seededPseudoRandomSample", (DL_FUNC) &seededPseudoRandomSample, 2},
  {NULL, NULL, 0}
};

/*
 * -- REGISTRATION OF THE .C ENTRY POINTS ---
 */
static const R_CMethodDef cMethods[] = {
  {"cleanupFxn", (DL_FUNC) &cleanupFxn, 0},
  {"shuffleRInt", (DL_FUNC) &shuffleRInt, 2},
  {"shuffleRRepl", (DL_FUNC) &shuffleRRepl, 2},
  {NULL, NULL, 0}
};

void R_init_SynExtend(DllInfo *info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
