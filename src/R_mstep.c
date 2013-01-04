/* This file contains a function R_mstep() called by R wraps m.step() in
   "R/fcn_m.step.r", and this function calls the relative functions mstep()
   in "src/emcluster.c" using in
   Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/10/19.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/emcluster.c".
   Input:
     **X: double[n, p], data matrix of n*p.
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     *pi: double[nclass], proportions of classes.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
     **Gamma: double[n, p], posterios matrix of n*p. 
   Output:
     *pi, **Mu, **LTSigma.
*/
void mstep(double **X,int n,int p,int k,double *pi,double **Mu,
           double **LTSigma,double **Gamma);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls mstep() in "src/emcluster.c" and is called by
   m.step() using .Call() in "R/fcn_m_step.r".
   Input:
     R_x: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.		# k
     R_Gamma: SEXP[R_n, R_p], posterios matrix of R_n*R_p. 
   Output:
     ret: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular sigma
                matrices.
*/
SEXP R_mstep(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_Gamma){
  /* Declare variables for calling C. */
  double **C_Gamma, **C_x, *C_pi, **C_Mu, **C_LTSigma;
  int *C_n, *C_p, *C_nclass;

  /* Declare variables for R's returning. */
  SEXP pi, Mu, LTSigma, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i, p_LTSigma;
  char *names[3] = {"pi", "Mu", "LTSigma"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  p_LTSigma = *C_p * (*C_p + 1) / 2;

  /* Allocate and protate storages. */
  PROTECT(pi = allocVector(REALSXP, *C_nclass));
  PROTECT(Mu = allocVector(REALSXP, *C_nclass * *C_p));
  PROTECT(LTSigma = allocVector(REALSXP, *C_nclass * p_LTSigma));
  PROTECT(ret = allocVector(VECSXP, 3));
  PROTECT(ret_names = allocVector(STRSXP, 3));

  i = 0;
  SET_VECTOR_ELT(ret, i++, pi);
  SET_VECTOR_ELT(ret, i++, Mu);
  SET_VECTOR_ELT(ret, i++, LTSigma);
  
  for(i = 0; i < 3; i++){
    SET_STRING_ELT(ret_names, i, mkChar(names[i])); 
  }
  setAttrib(ret, R_NamesSymbol, ret_names);

  /* Assign data. */
  C_Gamma = allocate_double_array(*C_n);
  C_x = allocate_double_array(*C_n);
  C_Mu = allocate_double_array(*C_nclass);
  C_LTSigma = allocate_double_array(*C_nclass);

  tmp_1 = REAL(R_Gamma);
  tmp_2 = REAL(R_x);
  for(i = 0; i < *C_n; i++){
    C_Gamma[i] = tmp_1;
    C_x[i] = tmp_2;
    tmp_1 += *C_nclass;
    tmp_2 += *C_p;
  }

  tmp_1 = REAL(Mu);
  tmp_2 = REAL(LTSigma);
  for(i = 0; i < *C_nclass; i++){
    C_Mu[i] = tmp_1;
    C_LTSigma[i] = tmp_2;
    tmp_1 += *C_p;
    tmp_2 += p_LTSigma;
  }

  C_pi = REAL(pi);

  /* Compute. */
  mstep(C_x, *C_n, *C_p, *C_nclass, C_pi, C_Mu, C_LTSigma, C_Gamma);

  /* Free memory and release protectation. */
  free(C_Gamma);
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(5);

  return(ret);
} /* End of R_mstep(). */

