/* This file contains several functions for counting iterations and time.
   Created: Wei-Chen Chen on 2009/06/05.
*/

#include <time.h>
#include "ss_mb_tool.h"

typedef struct emptr{
    double **C_X, *C_pi, **C_Mu, **C_LTSigma, *C_llhdval;
    int *C_n, *C_p, *C_nclass, *C_nc, *C_class;
    int *C_short_iter, *C_fixed_iter, *C_n_candidate, *C_EM_iter, *C_conv_iter;
    double *C_short_eps, *C_EM_eps, *C_conv_eps;
    int *C_lab, *C_labK;
    int *C_init_method;
    int C_protect_length;
} *EMPTR;

/* Functions in "it_R_emptr.c". */
SEXP create_emptr(SEXP R_X, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_short_iter, SEXP R_short_eps, SEXP R_fixed_iter, SEXP R_n_candidate,
    SEXP R_EM_iter, SEXP R_EM_eps, SEXP R_lab, SEXP R_labK,
    SEXP R_init_method,
    EMPTR emptr);

/* Functions in "it_init_EM.c". */
void init_EM(double **x, double *pi, double **Mu, double **LTSigma,
    double *llhdval,
    int n, int p, int nclass, int *nc, int *class,
    int short_iter, double short_eps, int fixed_iter, int n_candidate,
    int EM_iter, double EM_eps, int *conv_iter, double *conv_eps,
    int *lab, int labK,
    int init_method);

/* Functions in "it_emcluster.c". */
//void it_emcluster(int n, int p, int k, double *pi, double **X, double **Mu, 
//    double **LTSigma, int maxiter, double eps, double *llhdval);

/* Functions in "it_ss_emculster.c". */
//void it_ss_emcluster(int n, int p, int k, double *pi, double **X, double **Mu, 
//    double **LTSigma, int maxiter, double eps, double *llhdval,
//    int *lab);

