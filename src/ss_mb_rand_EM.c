/* This file is modified from "mb_rand_EM.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"

/* Modified shortems() for model-based initializer. */
void ss_mod_shortems_mb(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
    int *lab, int labK){
  int i, iter, totiter = 0, n_par = p * (p + 1) / 2;
  double *oldpi, **oldMu, **oldLTSigma, oldllh = -Inf, llhval;

  MAKE_VECTOR(oldpi, nclass);
  MAKE_MATRIX(oldMu, nclass, p);
  MAKE_MATRIX(oldLTSigma, nclass, n_par);

  do{
    ss_mb_init(X, n, p, nclass, oldpi, oldMu, oldLTSigma, lab, labK);

    iter = maxshortiter - totiter;
    if(fixed_iter > iter) fixed_iter = iter;
    iter = ss_mod_shortemcluster(n, p, nclass, oldpi, X, oldMu, oldLTSigma,
                                 fixed_iter, &llhval, lab);
    
    if(llhval >= oldllh){
      oldllh = llhval;
      cpy(oldMu, nclass, p, Mu);
      cpy(oldLTSigma, nclass, n_par,LTSigma);
      for(i = 0; i < nclass; i++) pi[i] = oldpi[i];
    }

    totiter += iter;
  } while(totiter < maxshortiter);

  FREE_MATRIX(oldMu);
  FREE_MATRIX(oldLTSigma);
  FREE_VECTOR(oldpi);
} /* End of ss_mod_shortems_mb(). */


/* This is the model-based rand.EM.
   This function is equal to mb_em_EM() with shorteps = Inf if fixed_iter = 1.
*/
int ss_mb_rand_EM(double **x, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma, double *llhdval, int shortiter, int fixed_iter,
    int *lab, int labK){
  if(nclass == 1){
    meandispersion_MLE(x, n, p, Mu[0], LTSigma[0]);
    *llhdval = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
               0.5 * n * p * log(2 * PI);
  } else {
    ss_mod_shortems_mb(n, p, nclass, pi, x, Mu, LTSigma, shortiter, fixed_iter,
                       lab, labK);
    ss_emcluster(n, p, nclass, pi, x, Mu, LTSigma, 1000, 0.0001, llhdval, lab);
  } 

  return 0;
} /* End of ss_mb_rand_EM(). */

