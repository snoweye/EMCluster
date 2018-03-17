/* This file contains a function R_emcluster() called by R wraps emcluster() in
   "R/fcn_emcluster.r", and this function calls the relative functions
   M_emcluster() in "src/M_emcluster.c" using in
   Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/10/12.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/ss_emcluster.c".
   Input:
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     *pi: double[nclass], proportions of classes.
     **X: double[n, p], data matrix of n*p.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
     maxiter: int[1], max iterations for emclust(), 1000 by default.	#em_iter
     eps: double[1], tolerance for emclust(), 1e-4 by default.		#em_eps
     *llhdval: double[1], log likelihood value.
     *lab: int[n], -1 for points with unknown clusters; 0,..,(labK-1) for known.
   Output:
     *pi, **Mu, **LTSigma, *llhdval, *conv_iter, *conv_eps.
*/
void ss_emcluster(int n,int p,int k,double *pi,double **X,double **Mu, 
    double **LTSigma,int maxiter,double eps,double *llhdval,
    int *conv_iter,double *conv_eps,int *lab);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls ss_emcluster() in "src/ss_emcluster.c" and is called by
   emcluster() using .Call() in "R/fcn_emcluster.r".
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
     R_em_iter: SEXP[1], max iterations for emclust(), 1000 by default.
     R_em_eps: SEXP[1], tolerance for emclust(), 1e-4 by default.
     R_lab: SEXP[1], -1 for points with unknown clusters;
                     0, .., (labK-1) for known.
   Output:
     ret: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular sigma
                matrices.
       llhdval: SEXP[1], log likelihood value.
       conv_iter: SEXP[1], convergent iterations.
       conv_eps: SEXP[1], convergent tolerance.
*/
SEXP ss_R_emcluster(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_p_LTSigma, SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_em_iter,
    SEXP R_em_eps, SEXP R_lab){
  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma, *C_llhdval, *C_em_eps, *C_conv_eps;
  int *C_n, *C_p, *C_nclass, *C_p_LTSigma, *C_em_iter, *C_conv_iter;
  int *C_lab;

  /* Declare variables for R's returning. */
  SEXP pi, Mu, LTSigma, llhdval, conv_iter, conv_eps, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i, j, tl;
  char *names[6] = {"pi", "Mu", "LTSigma", "llhdval", "conv.iter", "conv.eps"};

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
  PROTECT(conv_iter  = allocVector(INTSXP, 1));
  PROTECT(conv_eps  = allocVector(REALSXP, 1));
  PROTECT(ret = allocVector(VECSXP, 6));
  PROTECT(ret_names = allocVector(STRSXP, 6));

  i = 0;
  SET_VECTOR_ELT(ret, i++, pi);
  SET_VECTOR_ELT(ret, i++, Mu);
  SET_VECTOR_ELT(ret, i++, LTSigma);
  SET_VECTOR_ELT(ret, i++, llhdval);
  SET_VECTOR_ELT(ret, i++, conv_iter);
  SET_VECTOR_ELT(ret, i++, conv_eps);

  for(i = 0; i < 6; i++){
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
  C_em_iter = INTEGER(R_em_iter);
  C_em_eps = REAL(R_em_eps);
  C_lab = INTEGER(R_lab);
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
  ss_emcluster(*C_n, *C_p, *C_nclass, C_pi, C_x, C_Mu, C_LTSigma,
               *C_em_iter, *C_em_eps, C_llhdval,
               C_conv_iter, C_conv_eps, C_lab);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(8);

  return(ret);
} /* End of ss_R_emcluster(). */

