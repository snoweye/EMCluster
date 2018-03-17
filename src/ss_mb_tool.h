/* This file contains several functions to perform model-based initializers.
   Created: Wei-Chen Chen on 2009/03/04.
*/

#include "mb_tool.h"

/*-----------------------------------------------------------------------------
  The following functions are added by Wei-Chen Chen for EMCluster.
*/

/* Functions in "ss_M_emculster.c". */
void ss_norm_gamma(int n, int k, double **Gamma, double *pi, int *lab);
void ss_emcluster(int n, int p, int k, double *pi, double **X, double **Mu, 
    double **LTSigma, int maxiter, double eps, double *llhdval,
    int *conv_iter, double *conv_eps, int *lab);

/* Functions in "ss_emculster.c". */
void ss_estep(int n, int p, int k, double **X, double **Gamma, double *pi,
    double **Mu, double **LTSigma, int *lab);
void ss_emcluster_org(int n, int p, int k, double *pi, double **X, double **Mu, 
    double **LTSigma, int maxiter, double eps, double *llhdval, int *lab);

/* Functions in "ss_init_other.c". */
void ss_randomEMinit(double **x, int n, int p, int nclass, double *pi,
   double **Mu, double **LTSigma, int *lab, int labK, int nonlab_total,
   int *lab_index);
int ss_shortemcluster_org(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int maxiter, double eps, double *llhdval,
    int *lab);
int ss_shortemcluster(int n, int p, int k, double *pi, double **X, double **Mu,
    double **LTSigma, int maxiter, double eps, double *llhdval,
    int *conv_iter, double *conv_eps, int *lab);
int ss_shortems(int n, int p, int nclass, double *pi, double **X, double **Mu,  
    double **LTSigma, int maxshortiter, double shorteps,
    int *conv_iter, double *conv_eps, int *lab, int labK);

/* Functions in "ss_initials.c". */
void ss_assign(int n, int p, int k, double **X, double *pi, double **Mu,
    double **LTSigma, int *class, int *nc, int *lab);

/* Functions in "labInit.c". */
void labInitMus(int n, int p, int labK, double **x, int *lab, double **Mu); 


/* Functions in "ss_mb_em_EM.c". */
void ss_shortems_mb(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, double shorteps,
    int *conv_iter, double *conv_eps, int *lab, int labK);

/* Functions in "ss_mb_init.c". */
void ss_mb_init(double **X, int n, int p, int k, double *pi, double **Mu,
    double **LTSigma, int *lab, int labK);

/* Functions in "ss_mb_rand_EM.c". */
void ss_mod_shortems_mb(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
    int *conv_iter, double *conv_eps, int *lab, int labK);

/* Functions in "ss_mb_randomEMinit.c". */
int ss_mb_assign_closest(double **X, int n, int p, int nclass, double **Mu,
    double **LTSigma, int* clas, int *lab);
void ss_mb_randomEMinit(double **x, int n, int p, int nclass, double *pi,
    double **Mu, double **LTSigma, int *lab, int labK, int nonlab_total, int *lab_index);

/* Functions in "ss_rand_EM.c". */
int ss_mod_shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int fixed_iter, double *llhdval,
    int *conv_iter, double *conv_eps, int *lab);
void ss_mod_shortems(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
    int *conv_iter, double *conv_eps, int *lab, int labK);

/* Functions in "ss_ac_EM.c". */
// int ss_shortems_ac(int n, int p, int nclass, double *pi, double **X,
//     double **Mu, double **LTSigma, int maxshortiter, double shorteps,
//     int n_candidate, int *lab, int labK);
// void ss_mod_shortems_ac(int n, int p, int nclass, double *pi, double **X,
//     double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
//     int n_candidate, int *lab, int labK);

