/* This file contains a function R_M_emgroup() called by R wraps emgroup() in
   "R/fcn_emgroup.r", and this function calls the relative functions
   M_emgroup() in "src/M_emgroup.c" using in
   Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/10/07.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/M_emgroup.c".
   Input:
     **x: double[n, p], data matrix of n*p.
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     nclass: int[1], number of classes.
     *pi: double[nclass], proportions of classes.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
     *llhdval: double[1], log likelihood value.
     *nc: int[nclass], number of observations in each class.
     *class: int[n], class id's for all observations
             starting from 0 to (nclass - 1).
     alpha: double[1], 0.99 by default.
     em_iter: int[1], max iterations for emclust(), 1000 by default.
     em_eps: double[1], tolerance for emclust(), 1e-4 by default.
   Output:
     flag: int[1], returned value from M_emgroup().
     *pi, **Mu, **LTSigma, *llhdval, *nc, *class, *conv_iter, *conv_eps.
*/
int M_emgroup(double **x, int n, int p, int nclass, double *pi, double **Mu,
              double **LTSigma, double *llhdval,int *nc, int *class,
              double alpha, int em_iter, double em_eps,
              int *conv_iter, double *conv_eps);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls M_emgroup() in "src/M_emgroup.c" and is called by
   emgroup() using .Call() in "R/fcn_emgroup.r".
   Input:
     R_x: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.
     R_alpha: SEXP[1], 0.99 by default.
     R_em_iter: SEXP[1], max iterations for emclust(), 1000 by default.
     R_em_eps: SEXP[1], tolerance for emclust(), 1e-4 by default.
   Output:
     ret: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular sigma
                matrices.
       llhdval: SEXP[1], log likelihood value.
       nc: SEXP[R_nclass], number of observations in each class.
       class: SEXP[R_n], class id's for all observations
              starting from 0 to (R_nclass - 1).
       flag: SEXP[1], a returned value from M_emgroup() in "src/M_emgroup.c".
       conv_iter: SEXP[1], convergent iterations.
       conv_eps: SEXP[1], convergent tolerance.
*/
SEXP R_M_emgroup(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_alpha, SEXP R_em_iter, SEXP R_em_eps){
  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma, *C_llhdval, *C_alpha, *C_em_eps,
         *C_conv_eps;
  int *C_n, *C_p, *C_nclass, *C_nc, *C_class, *C_flag, *C_em_iter,
      *C_conv_iter;

  /* Declare variables for R's returning. */
  SEXP pi, Mu, LTSigma, llhdval, nc, class, flag, conv_iter, conv_eps,
       ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i, p_LTSigma;
  char *names[9] = {"pi", "Mu", "LTSigma", "llhdval", "nc", "class", "flag",
                    "conv.iter", "conv.eps"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  p_LTSigma = *C_p * (*C_p + 1) / 2;

  /* Allocate and protate storages. */
  PROTECT(pi = allocVector(REALSXP, *C_nclass));
  PROTECT(Mu = allocVector(REALSXP, *C_nclass * *C_p));
  PROTECT(LTSigma = allocVector(REALSXP, *C_nclass * p_LTSigma));
  PROTECT(llhdval = allocVector(REALSXP, 1));
  PROTECT(nc = allocVector(INTSXP, *C_nclass));
  PROTECT(class = allocVector(INTSXP, *C_n));
  PROTECT(flag = allocVector(INTSXP, 1));
  PROTECT(conv_iter = allocVector(INTSXP, 1));
  PROTECT(conv_eps = allocVector(REALSXP, 1));
  PROTECT(ret = allocVector(VECSXP, 9));
  PROTECT(ret_names = allocVector(STRSXP, 9));

  i = 0;
  SET_VECTOR_ELT(ret, i++, pi);
  SET_VECTOR_ELT(ret, i++, Mu);
  SET_VECTOR_ELT(ret, i++, LTSigma);
  SET_VECTOR_ELT(ret, i++, llhdval);
  SET_VECTOR_ELT(ret, i++, nc);
  SET_VECTOR_ELT(ret, i++, class);
  SET_VECTOR_ELT(ret, i++, flag);
  SET_VECTOR_ELT(ret, i++, conv_iter);
  SET_VECTOR_ELT(ret, i++, conv_eps);

  for(i = 0; i < 9; i++){
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
    tmp_2 += p_LTSigma;
  }

  C_pi = REAL(pi);
  C_llhdval = REAL(llhdval);
  C_nc = INTEGER(nc);
  C_class = INTEGER(class);
  C_flag = INTEGER(flag);
  C_alpha = REAL(R_alpha);
  C_em_iter = INTEGER(R_em_iter);
  C_em_eps = REAL(R_em_eps);
  C_conv_iter = INTEGER(conv_iter);
  C_conv_eps = REAL(conv_eps);

  /* Compute. */
  *C_flag = M_emgroup(C_x, *C_n, *C_p, *C_nclass, C_pi, C_Mu, C_LTSigma,
                      C_llhdval, C_nc, C_class,
                      *C_alpha, *C_em_iter, *C_em_eps,
                      C_conv_iter, C_conv_eps);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(11);

  return(ret);
} /* End of R_emgroup(). */

