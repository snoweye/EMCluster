/* This file contains functions called by R wraps in "R/fcn_init_other_shortem.r",
   and these functions call the relative functions in "src/init_other.c"
   using in Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2009/01/19.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/init_other.c".
   Input:
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.
     *pi: double[nclus], proportions of clusters.
     **X: double[n, p], data matrix of n*p.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclus, p * (p + 1) / 2], lower triangular sigma matrices.
     maxiter: int[1], number of short iterations, 500 by default.
     eps: double[1], epsilon for short iterations, 1e-2 by default.
   Output:
     *pi, **Mu, **LTSigma, *llhdval, iter, *conv_iter, *conv_eps.
*/
int shortemcluster(int n,int p,int k,double *pi,double **X,double **Mu,  
		   double **LTSigma,int maxiter,double eps,double *llhdval,
		   int *conv_iter,double *conv_eps);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls shortemcluster() in "src/init_other.c" and is called by
   shortemcluster() using .Call() in "R/fcn_init_other_shortem.r".
   Input:
     R_x: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.
     R_p_LTSigma: SEXP[1], dimersion of LTSigma, p * (p + 1) / 2.
     R_pi: SEXP[R_nclass], proportions of classes.
     R_Mu: SEXP[R_nclass, R_p], means of MVNs.
     R_LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular sigma
                matrices.
     R_maxiter: SEXP[1], number of short iterations, 500 by default.
     R_eps: SEXP[1], epsilon for short iterations, 1e-2 by default.
   Output:
     ret: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[nclass, R_p * (R_p + 1) / 2], lower triangular
                sigma matrices.
       llhdval: SEXP[1], log likelihood value.
       iter: SEXP[1], iterations used in short em.
       conv_iter: SEXP[1], convergent iterations.
       conv_eps: SEXP[1], convergent tolerance.
*/
SEXP R_shortemcluster(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_p_LTSigma, SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma,
    SEXP R_maxiter, SEXP R_eps){
  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma, *C_llhdval, *C_eps, *C_conv_eps;
  int *C_n, *C_p, *C_nclass, *C_p_LTSigma, *C_maxiter, *C_iter, *C_conv_iter;

  /* Declare variables for R's returning. */
  SEXP pi, Mu, LTSigma, llhdval, iter, conv_iter, conv_eps, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i, j, tl;
  char *names[7] = {"pi", "Mu", "LTSigma", "llhdval", "iter",
                    "conv.iter", "conv.eps"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  C_p_LTSigma = INTEGER(R_p_LTSigma);

  /* Allocate and protate storages. */
  PROTECT(pi = allocVector(REALSXP, *C_nclass));
  PROTECT(Mu = allocVector(REALSXP, *C_nclass * *C_p));
  PROTECT(LTSigma = allocVector(REALSXP, *C_nclass * *C_p_LTSigma));
  PROTECT(llhdval = allocVector(REALSXP, 1));
  PROTECT(iter = allocVector(INTSXP, 1));
  PROTECT(conv_iter  = allocVector(INTSXP, 1));
  PROTECT(conv_eps  = allocVector(REALSXP, 1));
  PROTECT(ret = allocVector(VECSXP, 7));
  PROTECT(ret_names = allocVector(STRSXP, 7));

  i = 0;
  SET_VECTOR_ELT(ret, i++, pi);
  SET_VECTOR_ELT(ret, i++, Mu);
  SET_VECTOR_ELT(ret, i++, LTSigma);
  SET_VECTOR_ELT(ret, i++, llhdval);
  SET_VECTOR_ELT(ret, i++, iter);
  SET_VECTOR_ELT(ret, i++, conv_iter);
  SET_VECTOR_ELT(ret, i++, conv_eps);

  for(i = 0; i < 7; i++){
    SET_STRING_ELT(ret_names, i, mkChar(names[i])); 
  }
  setAttrib(ret, R_NamesSymbol, ret_names);

  /* Assign data. */
  C_x = allocate_double_array(*C_n);
  C_Mu = allocate_double_array(*C_nclass);
  C_LTSigma = allocate_double_array(*C_nclass);

  tmp_1 = REAL(R_x);
  for(i = 0; i < *C_n; i++){
    C_x[i] = tmp_1;
    tmp_1 += *C_p;
  }

  tmp_1 = REAL(Mu);
  tmp_2 = REAL(LTSigma);
  for(i = 0; i < *C_nclass; i++){
    C_Mu[i] = tmp_1;
    C_LTSigma[i] = tmp_2;
    tmp_1 += *C_p;
    tmp_2 += *C_p_LTSigma;
  }

  C_pi = REAL(pi);
  C_llhdval = REAL(llhdval);
  C_iter = INTEGER(iter);
  C_maxiter = INTEGER(R_maxiter);
  C_eps = REAL(R_eps);
  C_conv_iter = INTEGER(conv_iter);
  C_conv_eps = REAL(conv_eps);

  /* Copy R objects to input oebjects for C. */
  tmp_1 = REAL(R_pi);
  for(i = 0; i < *C_nclass; i++){
    C_pi[i] = *(tmp_1 + i);
  }
  tl = 0;
  tmp_1 = REAL(R_Mu);
  for(i = 0; i < *C_nclass; i++){
    for(j = 0; j < *C_p; j++){
      C_Mu[i][j] = *(tmp_1 + tl++);
    }
  }
  tl = 0;
  tmp_1 = REAL(R_LTSigma);
  for(i = 0; i < *C_nclass; i++){
    for(j = 0; j < *C_p_LTSigma; j++){
      C_LTSigma[i][j] = *(tmp_1 + tl++);
    }
  }

  /* Compute. */
  *C_iter = shortemcluster(*C_n, *C_p, *C_nclass, C_pi, C_x, C_Mu, C_LTSigma,
                           *C_maxiter, *C_eps, C_llhdval,
                           C_conv_iter, C_conv_eps);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(9);

  return(ret);
} /* End of R_shortemcluster(). */

