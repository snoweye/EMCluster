/* This file contains three functions.
   1. R_lnlikelihood() is called by R wraps logL() in "R/fcn_dlmvnorm.r", and
      this function calls the relative functions lnlikelihood() in
      "src/dlmvnorm.c".
   2. R_mixllhd() is called by R wraps dmixmvn() in "R/fcn_dlmvnorm.r", and this
      function calls the relative functions mixllhd() in "src/dlmvnorm.c".
   3. R_dlmvnorm() called by R wraps dlmvn() in "R/fcn_dlmvnorm.r", and this
      function calls the relative functions dlmvnorm() in "src/dlmvnorm.c".
   All functions are using in Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/10/28.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/dlmvnorm.c".
   Input:
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     *pi: double[nclass], proportions of classes.
     **X: double[n, p], data matrix of n*p.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
   Output:
     llhd: double[1], logL.
*/
double lnlikelihood(int n,int p,int k,double *pi,double **X,double **Mu,
                    double **LTSigma);


/* This function is coded in "src/dlmvnorm.c".
   Input:
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     *x: double[p], data vector of length p.
     *pi: double[nclass], proportions of classes.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
   Output:
     dmixmvn: double[1], density of mixed multivariate normal distribution.
*/
double mixllhd(int p,int k,double *x,double *pi,double **Mu,double **LTSigma);


/* This function is coded in "src/dlmvnorm.c".
   Input:
     *x: double[p], data vector of length p.
     p: int[1], number of dimersions.
     *mu: double[p], means of MVN.
     *LTsigma: double[p * (p + 1) / 2], lower triangular sigma matrix.
   Output:
     exponent: double[1], log density of multivariate normal distribution.
*/
double dlmvnorm(double *x, int p, double *mu, double *LTsigma);


/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);




/* This function calls lnlikelihood() in "src/dlmvnorm.c" and is called by
   logL() using .Call() in "R/fcn_dlmvnorm.r".
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
   Output:
     logL: SEXP[1], a log likelihood value. 
*/
SEXP R_lnlikelihood(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_p_LTSigma, SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma){
  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma, *C_llhdval;
  int *C_n, *C_p, *C_nclass, *C_p_LTSigma;

  /* Declare variables for R's returning. */
  SEXP llhdval;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i;

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  C_p_LTSigma = INTEGER(R_p_LTSigma);

  /* Allocate and protate storages. */
  PROTECT(llhdval = allocVector(REALSXP, 1));

  /* Assign data. */
  C_x = allocate_double_array(*C_n);
  C_Mu = allocate_double_array(*C_nclass);
  C_LTSigma = allocate_double_array(*C_nclass);

  tmp_1 = REAL(R_x);
  for(i = 0; i < *C_n; i++){
    C_x[i] = tmp_1;
    tmp_1 += *C_p;
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
  C_llhdval = REAL(llhdval);

  /* Compute. */
  *C_llhdval = lnlikelihood(*C_n, *C_p, *C_nclass, C_pi, C_x, C_Mu, C_LTSigma);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(1);

  return(llhdval);
} /* End of R_lnlikelihood(). */




/* This function calls mixllhd() in "src/dlmvnorm.c" and is called by
   dmixmvn() using .Call() in "R/fcn_dlmvnorm.r".
   Input:
     R_x: SEXP[R_p], data vector of length R_p.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.		# k
     R_p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
     R_pi: SEXP[R_nclass], proportions of classes.
     R_Mu: SEXP[R_nclass, R_p], means of MVNs.
     R_LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular sigma
                matrices.
   Output:
     dmixmvn: SEXP[1], density of mixed multivariate normal distribution.
*/
SEXP R_mixllhd(SEXP R_x, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma){
  /* Declare variables for calling C. */
  double *C_x, *C_pi, **C_Mu, **C_LTSigma, *C_dmixmvn;
  int *C_p, *C_nclass, *C_p_LTSigma;

  /* Declare variables for R's returning. */
  SEXP dmixmvn;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i;

  /* Set initial values. */
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  C_p_LTSigma = INTEGER(R_p_LTSigma);

  /* Allocate and protate storages. */
  PROTECT(dmixmvn = allocVector(REALSXP, 1));

  /* Assign data. */
  C_Mu = allocate_double_array(*C_nclass);
  C_LTSigma = allocate_double_array(*C_nclass);

  tmp_1 = REAL(R_Mu);
  tmp_2 = REAL(R_LTSigma);
  for(i = 0; i < *C_nclass; i++){
    C_Mu[i] = tmp_1;
    C_LTSigma[i] = tmp_2;
    tmp_1 += *C_p;
    tmp_2 += *C_p_LTSigma;
  }

  C_x = REAL(R_x);
  C_pi = REAL(R_pi);
  C_dmixmvn = REAL(dmixmvn);

  /* Compute. */
  *C_dmixmvn = mixllhd(*C_p, *C_nclass, C_x, C_pi, C_Mu, C_LTSigma);

  /* Free memory and release protectation. */
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(1);

  return(dmixmvn);
} /* End of R_mixllhd(). */




/* This function calls dlmvnorm() in "src/dlmvnorm.c" and is called by
   dlmvn() using .Call() in "R/fcn_dlmvnorm.r".
   Input:
     R_x: SEXP[R_p], data vector of length R_p.
     R_p: SEXP[1], number of dimersions.
     R_p_LTsigma: SEXP[1], length of LTsigma, p * (p + 1) / 2.
     R_mu: SEXP[R_p], means of MVN.
     R_LTsigma: SEXP[R_p * (R_p + 1) / 2], lower triangular sigma matrix.
   Output:
     dlmvn: SEXP[1], log density of multivariate normal distribution.
*/
SEXP R_dlmvnorm(SEXP R_x, SEXP R_p, SEXP R_p_LTsigma, SEXP R_mu,
    SEXP R_LTsigma){
  /* Declare variables for calling C. */
  double *C_x, *C_mu, *C_LTsigma, *C_dlmvn;
  int *C_p, *C_p_LTsigma;

  /* Declare variables for R's returning. */
  SEXP dlmvn;

  /* Set initial values. */
  C_x = REAL(R_x);
  C_p = INTEGER(R_p);
  C_p_LTsigma = INTEGER(R_p_LTsigma);
  C_mu = REAL(R_mu);
  C_LTsigma = REAL(R_LTsigma);

  /* Allocate and protate storages. */
  PROTECT(dlmvn = allocVector(REALSXP, 1));

  /* Assign data. */
  C_dlmvn = REAL(dlmvn);

  /* Compute. */
  *C_dlmvn = dlmvnorm(C_x, *C_p, C_mu, C_LTsigma);

  /* Free memory and release protectation. */
  UNPROTECT(1);

  return(dlmvn);
} /* End of R_dlmvnorm(). */

