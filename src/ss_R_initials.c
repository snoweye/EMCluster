/* This file is modified from "R_initials.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/ss_initials.c".
   Input:
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     **X: double[n, p], data matrix of n*p.
     *pi: double[nclass], proportions of classes.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
     *lab: int[n], -1 for points with unknown clusters; 0,..,(labK-1) for known.
   Output:
     *class: int[n], class id's for all observations
             starting from 0 to (nclass - 1).
     *nc: int[nclass], number of observations in each class.
*/
void ss_assign(int n, int p, int k, double **X, double *pi, double **Mu,
    double **LTSigma, int *class, int *nc, int *lab);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls ss_assign() in "src/ss_initials.c" and is called by
   assign.class() using .Call() in "R/fcn_initials.r".
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
     R_lab: SEXP[1], -1 for points with unknown clusters;
                     0, .., (labK-1) for known.
   Output:
     ret: a list contains
       nc: SEXP[R_nclass], number of observations in each class.
       class: SEXP[R_n], class id's for all observations
              starting from 0 to (R_nclass - 1).
*/
SEXP ss_R_assign(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_lab){
  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma;
  int *C_n, *C_p, *C_nclass, *C_p_LTSigma, *C_nc, *C_class;
  int *C_lab;

  /* Declare variables for R's returning. */
  SEXP nc, class, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1, *tmp_2;
  int i;
  char *names[2] = {"nc", "class"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_nclass = INTEGER(R_nclass);
  C_p_LTSigma = INTEGER(R_p_LTSigma);

  /* Allocate and protate storages. */
  PROTECT(nc = allocVector(INTSXP, *C_nclass));
  PROTECT(class = allocVector(INTSXP, *C_n));
  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(ret_names = allocVector(STRSXP, 2));

  i = 0;
  SET_VECTOR_ELT(ret, i++, nc);
  SET_VECTOR_ELT(ret, i++, class);

  for(i = 0; i < 2; i++){
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

  tmp_1 = REAL(R_Mu);
  tmp_2 = REAL(R_LTSigma);
  for(i = 0; i < *C_nclass; i++){
    C_Mu[i] = tmp_1;
    C_LTSigma[i] = tmp_2;
    tmp_1 += *C_p;
    tmp_2 += *C_p_LTSigma;
  }

  C_pi = REAL(R_pi);
  C_nc = INTEGER(nc);
  C_class = INTEGER(class);
  C_lab = INTEGER(R_lab);

  /* Compute. */
  ss_assign(*C_n, *C_p, *C_nclass, C_x, C_pi, C_Mu, C_LTSigma, C_class, C_nc,
            C_lab);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(4);

  return(ret);
} /* End of ss_R_assign(). */

