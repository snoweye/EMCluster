/* This file contains two functions R_assign() R_meandispersion()
   called by R wraps assign.class() and meandispersion()
   in "R/fcn_initials.r", and these functions call the
   relative functions assign() and meandispersion()
   in "src/initials.c" using in Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2008/12/05.
*/

#include<R.h>
#include<Rinternals.h>


/* This function is coded in "src/initials.c".
   Input:
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     k: int[1], number of classes.		# nclass
     **X: double[n, p], data matrix of n*p.
     *pi: double[nclass], proportions of classes.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclass, p * (p + 1) / 2], lower triangular sigma
                matrices.
   Output:
     *class: int[n], class id's for all observations
             starting from 0 to (nclass - 1).
     *nc: int[nclass], number of observations in each class.
*/
void assign(int n, int p,int k,double **X,double *pi,double **Mu,
	   double **LTSigma,int *class,int *nc);

/* This function is coded in "src/initials.c".
   Input:
     **x: double[n, p], data matrix of n*p.
     n: int[1], number of observations.
     p: int[1], number of dimersions.
   Output:
     *mu: double[p], means of MVNs.
     *ltsigma: double[p * (p + 1) / 2], lower triangular sigma matrices.
*/
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);
void meandispersion_MLE(double **x, int n, int p, double *mu, double *ltsigma);
void meandispersion_MME(double **x, int n, int p, double *mu, double *ltsigma);

/* Allocate a pointer array with double precision. See "src/R_tool.c". */
double** allocate_double_array(int n);


/* This function calls assign() in "src/initials.c" and is called by
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
   Output:
     ret: a list contains
       nc: SEXP[R_nclass], number of observations in each class.
       class: SEXP[R_n], class id's for all observations
              starting from 0 to (R_nclass - 1).
*/
SEXP R_assign(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma){
  /* Declare variables for calling C. */
  double **C_x, *C_pi, **C_Mu, **C_LTSigma;
  int *C_n, *C_p, *C_nclass, *C_p_LTSigma, *C_nc, *C_class;

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

  /* Compute. */
  assign(*C_n, *C_p, *C_nclass, C_x, C_pi, C_Mu, C_LTSigma, C_class, C_nc);

  /* Free memory and release protectation. */
  free(C_x);
  free(C_Mu);
  free(C_LTSigma);
  UNPROTECT(4);

  return(ret);
} /* End of R_assign(). */




/* This function calls meandispersion() in "src/initials.c" and is called by
   meandispersion() using .Call() in "R/fcn_initials.r".
   Input:
     R_x: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_type: SEXP[1], 0 for original version, 1 for MLE, 2 for MME.
   Output:
     ret: a list contains
       mu: SEXP[R_p], means of MVNs.
       ltsigma: SEXP[R_p * (R_p + 1) / 2], lower triangular sigma matrices.
*/
SEXP R_meandispersion(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_type){
  /* Declare variables for calling C. */
  double **C_x, *C_mu, *C_ltsigma;
  int *C_n, *C_p, *C_type;

  /* Declare variables for R's returning. */
  SEXP mu, ltsigma, ret, ret_names;

  /* Declare variables for processing. */
  double *tmp_1;
  int i;
  char *names[2] = {"mu", "ltsigma"};

  /* Set initial values. */
  C_n = INTEGER(R_n);
  C_p = INTEGER(R_p);
  C_type = INTEGER(R_type);

  /* Allocate and protate storages. */
  PROTECT(mu = allocVector(REALSXP, *C_p));
  PROTECT(ltsigma = allocVector(REALSXP, *C_p * (*C_p + 1) / 2));
  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(ret_names = allocVector(STRSXP, 2));

  i = 0;
  SET_VECTOR_ELT(ret, i++, mu);
  SET_VECTOR_ELT(ret, i++, ltsigma);

  for(i = 0; i < 2; i++){
    SET_STRING_ELT(ret_names, i, mkChar(names[i])); 
  }
  setAttrib(ret, R_NamesSymbol, ret_names);

  /* Assign data. */
  C_x = allocate_double_array(*C_n);

  tmp_1 = REAL(R_x);
  for(i = 0; i < *C_n; i++){
    C_x[i] = tmp_1;
    tmp_1 += *C_p;
  }

  C_mu = REAL(mu);
  C_ltsigma = REAL(ltsigma);

  /* Compute. */
  switch(*C_type){
    case 1:
      meandispersion_MLE(C_x, *C_n, *C_p, C_mu, C_ltsigma);
      break;
    case 2:
      meandispersion_MME(C_x, *C_n, *C_p, C_mu, C_ltsigma);
      break;
    default:
      meandispersion(C_x, *C_n, *C_p, C_mu, C_ltsigma);
  }

  /* Free memory and release protectation. */
  free(C_x);
  UNPROTECT(4);

  return(ret);
} /* End of R_meandispersion(). */

