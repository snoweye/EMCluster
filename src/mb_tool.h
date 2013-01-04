/* This file contains several functions to perform model-based initializers.
   Created: Wei-Chen Chen on 2009/03/04.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>


/* These are all required for EMCluster and defined at other files
   provided by Dr. Ranjan Maitra and Volodymyr Melnykov.
*/
#include "array.h"
#include "mat_vec.h"
#include "order.h"
#include "quantile.h"
#define Inf 1e+140

int srswor(int n, int k, int *ranordr);
int classify(double *X,int p,int k,double *pi, double **Mu, double **LTSigma);
int assign_closest(double *X, int p, int nclass, double **Mu);
double findzero(double ax, double bx, double (*f)(double x, void *info),
    void *info, double *Tol, int *Maxit);

double determinant(double *LTSigma,int n);
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);
double dlmvnorm(double *x, int p, double *mu, double *LTsigma);
double lnlikelihood(int n,int p,int k,double *pi,double **X,double **Mu,
    double **LTSigma);
double mixllhd(int p,int k,double *x,double *pi,double **Mu,double **LTSigma);
int emcluster_org(int n,int p,int k,double *pi,double **X,double **Mu, 
    double **LTSigma,int maxiter,double eps,double *llhdval);
void estep(int n,int p,int k,double **X,double **Gamma,double *pi,double **Mu, 
    double **LTSigma);
void mstep(double **X,int n,int p,int k,double *pi,double **Mu,
    double **LTSigma,double **Gamma);
void assign(int n, int p,int k,double **X,double *pi,double **Mu,
    double **LTSigma,int *class,int *nc);
int shortemcluster_org(int n,int p,int k,double *pi,double **X,double **Mu,  
    double **LTSigma,int maxiter,double eps,double *llhdval);
int shortems(int n,int p,int nclass,double *pi,double **X,double **Mu,  
    double **LTSigma,int maxshortiter,double shorteps);

int initials(double **x,int n,int p,int nclass,int *nc,
    double **Mu,double **LTSigma,int *class);
int randomEMinit(double **x,int n,int p,int nclass,double *pi,
    double **Mu,double **LTSigma);
int starts_via_svd(int n,int m,double **Mu,double **x,int nclus,int *ningrp,
    double *pi,int *grpids,double **LTSigma,double alpha, int llhdnotW);


/*-----------------------------------------------------------------------------
  The following functions are added by Wei-Chen Chen for EMCluster.
*/

/* Function in "M_emculster". */
double lnlikelihood_gamma(int n, int k, double **Gamma, double *pi);
void estep_gamma(int n, int p, int k, double **X, double **Gamma,
    double **Mu, double **LTSigma);
void estep_unnorm_dlmvn(int n, int p, int k, double **X, double **Gamma,
    double *pi, double **Mu, double **LTSigma);
void norm_gamma(int n, int k, double **Gamma, double *pi);
void emcluster(int n, int p, int k, double *pi, double **X, double **Mu, 
    double **LTSigma, int maxiter, double eps, double *llhdval);

/* Function in "M_init_other.c". */
int shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int maxiter, double eps, double *llhdval);

/* Functions in "meandispersion.c". */
void meandispersion_MLE(double **x, int n, int p, double *mu, double *ltsigma);
void meandispersion_MME(double **x, int n, int p, double *mu, double *ltsigma);
void est_ltsigma_mle_given_mu(double **x, int n, int p, double *mu,
    double *ltsigma);

/* Functions in "rand_EM.c". */
int mod_shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int fixed_iter, double *llhdval);
void mod_shortems(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter);

/* Functions in "Rtool.c". */
double** allocate_double_array(int n);

/* Functions in "M_emgroup.c". */
int M_emgroup(double **x,int n,int p,int nclass,double *pi,double **Mu,
              double **LTSigma,double *llhdval,int *nc,int *class,
              double alpha, int em_iter, double em_eps);

/* Functions in "mb_em_EM.c". */
void shortems_mb(int n, int p, int nclass, double *pi, double **X, double **Mu,  
    double **LTSigma, int maxshortiter, double shorteps);

/* Functions in "mb_init.c". */
void cut_sub(double **X, int n, int p, int G, int min_n, double lambda,
    double *prob, double **Mu, double **LTSigma);
void mb_init(double **X, int n, int p, int k, double *pi, double **Mu,
    double **LTSigma);

/* Functions in "mb_rand_EM.c". */
void mod_shortems_mb(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter);

/* Function in "mb_randomEMinit.c". */
int mb_assign_closest(double **X, int n, int p, int nclass, double **Mu,
    double **LTSigma, int* clas);
void mb_randomEMinit(double **x, int n, int p, int nclass, double *pi,
    double **Mu, double **LTSigma);

/* Functions in "mb_tool.c". */
typedef struct Param{
  double lower_bound, upper_bound, n_nclass, tol;
  int maxit;
} PARAM; /* End of Param. */
double mb_median(int n, double *x);
double mb_quantile(int n, double *x, double p);
double fcn(double lambda, void *pt_param);
double find_lambda(void *pt_param);

/* Functions in "ac_EM.c". */
int shortems_ac(int n, int p, int nclass, double *pi, double **X, double **Mu,  
    double **LTSigma, int maxshortiter, double shorteps, int n_candidate);
void mod_shortems_ac(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
    int n_candidate);

