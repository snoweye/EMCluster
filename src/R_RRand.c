/* This file contains a function R_RRand() called by R wraps rrand() in
   "R/fcn_rrand.r", and this function calls the relative functions RRand() in
   "src/RRand.c" using in
   Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/10/27.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/RRand.c".
   Input:
     *N: int[1], number of observations.
     *TRUK: int[1], number of true clusters.
     *PREDK: int[1], number of predicted clusters.
     *trcl: int[N], true cluster ids.
     *prcl: int[N], predicted cluster ids.
     *Rand: double[1], Rand index.
     *adjRand: double[1], adjust Rand index.
     *Eindex: double[1], Eindex.
   Output:
     *Rand, *adjRand, *Eindex.
*/
void RRand(int *N, int *TRUK, int *PREDK, int *trcl, int *prcl,
	   double *Rand, double *adjRand, double *Eindex);


/* This function calls RRand() in "src/RRand.c" and is called by
   rrand() using .Call() in "R/fcn_rrand.r".
   Input:
     R_N: SEXP[1], number of observations.
     R_TRUK: SEXP[1], number of true clusters.
     R_PREDK: SEXP[1], number of predicted clusters.
     R_trcl: SEXP[N], true cluster ids.
     R_prcl: SEXP[N], predicted cluster ids.
     R_Rand: SEXP[1], Rand index.
     R_adjRand: SEXP[1], adjust Rand index.
     R_Eindex: SEXP[1], Eindex.
   Output:
     ret: a list contains
       Rand: SEXP[1], Rand index.
       adjRand: SEXP[1], adjust Rand index.
       Eindex: SEXP[1], Eindex.
*/
SEXP R_RRand(SEXP R_N, SEXP R_TRUK, SEXP R_PREDK, SEXP R_trcl, SEXP R_prcl){
  /* Declare variables for calling C. */
  double *C_Rand, *C_adjRand, *C_Eindex;
  int *C_N, *C_TRUK, *C_PREDK, *C_trcl, *C_prcl;

  /* Declare variables for R's returning. */
  SEXP Rand, adjRand, Eindex, ret, ret_names;

  /* Declare variables for processing. */
  int i;
  char *names[3] = {"Rand", "adjRand", "Eindex"};
  
  /* Set initial values. */
  C_N = INTEGER(R_N);
  C_TRUK = INTEGER(R_TRUK);
  C_PREDK = INTEGER(R_PREDK);
  C_trcl = INTEGER(R_trcl);
  C_prcl = INTEGER(R_prcl);
  
  PROTECT(Rand = allocVector(REALSXP, 1));
  PROTECT(adjRand = allocVector(REALSXP, 1));
  PROTECT(Eindex = allocVector(REALSXP, 1));
  PROTECT(ret = allocVector(VECSXP, 3));
  PROTECT(ret_names = allocVector(STRSXP, 3));
  
  i = 0;
  SET_VECTOR_ELT(ret, i++, Rand);
  SET_VECTOR_ELT(ret, i++, adjRand);
  SET_VECTOR_ELT(ret, i++, Eindex);

  for(i = 0; i < 3; i++){
    SET_STRING_ELT(ret_names, i, mkChar(names[i])); 
  }
  setAttrib(ret, R_NamesSymbol, ret_names);

  /* Assign data. */
  C_Rand = REAL(Rand);
  C_adjRand = REAL(adjRand);
  C_Eindex = REAL(Eindex);

  /* Compute. */
  RRand(C_N, C_TRUK, C_PREDK, C_trcl, C_prcl, C_Rand, C_adjRand, C_Eindex);

  /* Free memory and release protectation. */
  UNPROTECT(5);

  return(ret);
} /* End of R_RRand(). */

