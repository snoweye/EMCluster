/* This file contains a function to create emptr for R.

   Writen: Wei-Chen Chen on 2009/04/27.
*/

#include<R.h>
#include<Rinternals.h>
#include "it_tool.h"

/*
   Input:
     R_X: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.
     C_*: points for other functions.
   Output:
     ret: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[nclass, R_p * (R_p + 1) / 2], lower triangular
                sigma matrices.
       llhdval: SEXP[1], log likelihood value.
       nc: SEXP[R_nclass], number of observations in each class.
       class: SEXP[R_n], class id's for all observations
              starting from 0 to (R_nclass - 1).
*/
SEXP create_emptr(SEXP R_X, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_short_iter, SEXP R_short_eps, SEXP R_fixed_iter, SEXP R_n_candidate,
    SEXP R_EM_iter, SEXP R_EM_eps, SEXP R_lab, SEXP R_labK,
    SEXP R_init_method,
    EMPTR emptr){
  /* Declare variables for R's returning. */
  SEXP pi, Mu, LTSigma, llhdval, nc, class,
       ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i, p_LTSigma;
  int ret_length = 6;
  char *names[6] = {"pi", "Mu", "LTSigma", "llhdval", "nc", "class"};

  emptr->C_protect_length = ret_length + 2;

  /* Set initial values. */
  emptr->C_n = INTEGER(R_n);
  emptr->C_p = INTEGER(R_p);
  emptr->C_nclass = INTEGER(R_nclass);
  p_LTSigma = *emptr->C_p * (*emptr->C_p + 1) / 2;

  /* Allocate and protate storages. */
  PROTECT(pi = allocVector(REALSXP, *emptr->C_nclass));
  PROTECT(Mu = allocVector(REALSXP, *emptr->C_nclass * *emptr->C_p));
  PROTECT(LTSigma = allocVector(REALSXP, *emptr->C_nclass * p_LTSigma));
  PROTECT(llhdval = allocVector(REALSXP, 1));
  PROTECT(nc = allocVector(INTSXP, *emptr->C_nclass));
  PROTECT(class = allocVector(INTSXP, *emptr->C_n));
  PROTECT(ret = allocVector(VECSXP, ret_length));
  PROTECT(ret_names = allocVector(STRSXP, ret_length));

  i = 0;
  SET_VECTOR_ELT(ret, i++, pi);
  SET_VECTOR_ELT(ret, i++, Mu);
  SET_VECTOR_ELT(ret, i++, LTSigma);
  SET_VECTOR_ELT(ret, i++, llhdval);
  SET_VECTOR_ELT(ret, i++, nc);
  SET_VECTOR_ELT(ret, i++, class);

  for(i = 0; i < ret_length; i++){
    SET_STRING_ELT(ret_names, i, mkChar(names[i])); 
  }
  setAttrib(ret, R_NamesSymbol, ret_names);

  /* Assign data. */
  emptr->C_X = allocate_double_array(*emptr->C_n);
  emptr->C_Mu = allocate_double_array(*emptr->C_nclass);
  emptr->C_LTSigma = allocate_double_array(*emptr->C_nclass);

  tmp_1 = REAL(R_X);
  for(i = 0; i < *emptr->C_n; i++){
    emptr->C_X[i] = tmp_1;
    tmp_1 += *emptr->C_p;
  }
  tmp_1 = REAL(Mu);
  tmp_2 = REAL(LTSigma);
  for(i = 0; i < *emptr->C_nclass; i++){
    emptr->C_Mu[i] = tmp_1;
    emptr->C_LTSigma[i] = tmp_2;
    tmp_1 += *emptr->C_p;
    tmp_2 += p_LTSigma;
  }

  emptr->C_pi = REAL(pi);
  emptr->C_llhdval = REAL(llhdval);
  emptr->C_nc = INTEGER(nc);
  emptr->C_class = INTEGER(class);
  emptr->C_short_iter = INTEGER(R_short_iter);
  emptr->C_short_eps = REAL(R_short_eps);
  emptr->C_fixed_iter = INTEGER(R_fixed_iter);
  emptr->C_n_candidate = INTEGER(R_n_candidate);
  emptr->C_EM_iter = INTEGER(R_EM_iter);
  emptr->C_EM_eps = REAL(R_EM_eps);
  emptr->C_lab = INTEGER(R_lab);
  emptr->C_labK = INTEGER(R_labK);
  emptr->C_init_method = INTEGER(R_init_method);

  return(ret);
} /* End of create_emptr(). */

