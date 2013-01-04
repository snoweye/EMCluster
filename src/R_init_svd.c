/* This file contains a function R_starts_via_svd() called by R wraps
   starts.via.svd() in "R/fcn_init_svd.r", and this function calls the
   relative functions starts_via_svd() in "src/init_svd.c" using in
   Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/09/28.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/init_svd.c".
   Input:
     n: int[1], number of observations.
     m: int[1], number of dimersions.
     **Mu: double[nclass, p], means of MVNs.
     **x: double[n, m], data matrix of n*m.
     nclus: int[1], number of clusters.
     *ningrp: int[nclus], number of observations in each cluster.
     *pi: double[nclus], proportions of clusters.
     *grpids: int[n], class id's for all observations starting
              from 0 to (nclus - 1).
     **LTSigma: double[nclus, m * (m + 1) / 2], lower triangular sigma matrices.
     alpha: double[1], 0.99.
     llhdnotW: int[1], 0 for kmeans, 1 for EM.
   Output:
     ind: int[1], a returned index from starts_via_svd() and starts_svd().
     **Mu, *ningrp, *pi, *grpids, **LTSigma,.
*/
int starts_via_svd(int n, int m, double **Mu, double **x, int nclus,
                   int *ningrp, double *pi, int *grpids, double **LTSigma,
                   double alpha, int llhdnotW);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls starts_via_svd() in "src/init_svd.c" and is called by
   starts.via.svd() using .Call() in "R/fcn_init_svd.r".
   Input:
     R_x: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.		(m)
     R_nclass: SEXP[1], number of classes.		(nclus)
     R_method: SEXP[1], 0 for kmeans, 1 for EM.		(llhdnotW)
     R_alpha: SEXP[1], 0.99 by default.
   Output:
     ret: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[R_nclass, R_p * (R_p + 1) / 2], lower triangular
                sigma matrices.
       nc: SEXP[R_nclass], number of observations in each class.	(ningrp)
       class: SEXP[R_n], class id's for all observations
              starting from 0 to (R_nclass - 1).			(grpids)
       flag: SEXP[1], a returned value from starts_via_svd() in
             "src/init_svd.c".						(ind)
*/
SEXP R_starts_via_svd(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_method, SEXP R_alpha){

  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma, *C_alpha;
  int *C_n, *C_p, *C_nclass, *C_nc, *C_class, *C_method, *C_flag;

  /* Declare variables for R's returning. */
  SEXP pi, Mu, LTSigma, nc, class, flag, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i, p_LTSigma;
  char *names[6] = {"pi", "Mu", "LTSigma", "nc", "class","flag"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  C_method = INTEGER(R_method);
  C_alpha = REAL(R_alpha);
  p_LTSigma = *C_p * (*C_p + 1) / 2;

  /* Allocate and protate storages. */
  PROTECT(pi = allocVector(REALSXP, *C_nclass));
  PROTECT(Mu = allocVector(REALSXP, *C_nclass * *C_p));
  PROTECT(LTSigma = allocVector(REALSXP, *C_nclass * p_LTSigma));
  PROTECT(nc = allocVector(INTSXP, *C_nclass));
  PROTECT(class = allocVector(INTSXP, *C_n));
  PROTECT(flag = allocVector(INTSXP, 1));
  PROTECT(ret = allocVector(VECSXP, 6));
  PROTECT(ret_names = allocVector(STRSXP, 6));

  i = 0;
  SET_VECTOR_ELT(ret, i++, pi);
  SET_VECTOR_ELT(ret, i++, Mu);
  SET_VECTOR_ELT(ret, i++, LTSigma);
  SET_VECTOR_ELT(ret, i++, nc);
  SET_VECTOR_ELT(ret, i++, class);
  SET_VECTOR_ELT(ret, i++, flag);

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
    tmp_2 += p_LTSigma;
  }

  C_pi = REAL(pi);
  C_nc = INTEGER(nc);
  C_class = INTEGER(class);
  C_flag = INTEGER(flag);

  /* Compute. */
  *C_flag = starts_via_svd(*C_n, *C_p, C_Mu, C_x, *C_nclass, C_nc,
                           C_pi, C_class, C_LTSigma, *C_alpha, *C_method);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(8);

  return(ret);
} /* End of R_starts_via_svd(). */

