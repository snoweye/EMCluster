/* This file contains a function R_estep() called by R wraps e.step() in
   "R/fcn_e_step.r", and this function calls the relative functions estep() in
   "src/emcluster.c" using in
   Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/10/19.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/emcluster.c".
   Input:
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     **X: double[n, p], data matrix of n*p.
     **Gamma: double[n, p], posterios matrix of n*p. 
     *pi: double[nclass], proportions of classes.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
     norm: int[1], normalized.
   Output:
     **Gamma.
*/
void estep(int n, int p, int k, double **X, double **Gamma, double *pi,
    double **Mu, double **LTSigma);
void estep_unnorm_dlmvn(int n, int p, int k, double **X, double **Gamma,
    double *pi, double **Mu, double **LTSigma);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls estep() in "src/emcluster.c" and is called by
   e.step() using .Call() in "R/fcn_e_step.r".
   Input:
     R_x: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.		# k
     R_p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
     R_pi: SEXP[R_nclass], proportions of classes.
     R_Mu: SEXP[R_nclass, R_p], means of MVNs.
     R_LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular sigma
                matrices.
     R_norm: SEXP[1], normalized.
   Output:
     ret: a list contains
       Gamma: SEXP[R_n, R_p], posterios matrix of R_n*R_p. 
*/
SEXP R_estep(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_norm){
  /* Declare variables for calling C. */
  double **C_Gamma, **C_x, *C_pi, **C_Mu, **C_LTSigma;
  int *C_n, *C_p, *C_nclass, *C_p_LTSigma, *C_norm;

  /* Declare variables for R's returning. */
  SEXP Gamma, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i;
  char *names[1] = {"Gamma"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  C_p_LTSigma = INTEGER(R_p_LTSigma);

  /* Allocate and protate storages. */
  PROTECT(Gamma = allocVector(REALSXP, *C_n * *C_nclass));
  PROTECT(ret = allocVector(VECSXP, 1));
  PROTECT(ret_names = allocVector(STRSXP, 1));

  SET_VECTOR_ELT(ret, 0, Gamma);
  SET_STRING_ELT(ret_names, 0, mkChar(names[0])); 
  setAttrib(ret, R_NamesSymbol, ret_names);

  /* Assign data. */
  C_Gamma = allocate_double_array(*C_n);
  C_x = allocate_double_array(*C_n);
  C_Mu = allocate_double_array(*C_nclass);
  C_LTSigma = allocate_double_array(*C_nclass);

  tmp_1 = REAL(Gamma);
  tmp_2 = REAL(R_x);
  for(i = 0; i < *C_n; i++){
    C_Gamma[i] = tmp_1;
    C_x[i] = tmp_2;
    tmp_1 += *C_nclass;
    tmp_2 += *C_p;
  }

  tmp_1 = REAL(R_Mu);
  tmp_2 = REAL(R_LTSigma);
  for(i = 0; i < *C_nclass; i++){
    C_Mu[i] = tmp_1;
    C_LTSigma[i] = tmp_2;
    tmp_1 += *C_p;
    tmp_2 += *C_p_LTSigma;
  }

  C_pi = REAL(R_pi);
  C_norm = INTEGER(R_norm);

  /* Compute. */
  if(*C_norm == 1){
    estep(*C_n, *C_p, *C_nclass, C_x, C_Gamma, C_pi, C_Mu, C_LTSigma);
  } else{
    estep_unnorm_dlmvn(*C_n, *C_p, *C_nclass, C_x, C_Gamma, C_pi, C_Mu,
                       C_LTSigma);
  }

  /* Free memory and release protectation. */
  free(C_Gamma);
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(3);

  return(ret);
} /* End of R_estep(). */

