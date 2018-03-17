/* This file contains several functions modified from Dr. Maitra's original
   codes in "init_other.c" and "M_init_other.c" to perform the modified
   rand.EM method.
   Created: Wei-Chen Chen on 2009/02/07.
*/

#include "mb_tool.h"

/* Modified shortemcluster(). */
int mod_shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int fixed_iter, double *llhdval,
    int *conv_iter, double *conv_eps){
  int iter;
  double **gamm;

  MAKE_MATRIX(gamm,n,k);

  iter=0;
  do{
    estep(n, p, k, X, gamm, pi, Mu, LTSigma);
    mstep(X, n, p, k, pi, Mu, LTSigma, gamm);
    iter++;
  } while(iter < fixed_iter);
  *conv_iter = iter;
  *conv_eps = -1.0;

  *llhdval = lnlikelihood(n, p, k, pi, X, Mu, LTSigma);
  FREE_MATRIX(gamm);

  return iter;
} /* End of mod_shortemcluster(). */


/* Modified shortems(). */
void mod_shortems(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
    int *conv_iter, double *conv_eps){
  int i, iter, totiter = 0, n_par = p * (p + 1) / 2;
  double *oldpi, **oldMu, **oldLTSigma, oldllh = -Inf, llhval;

  MAKE_VECTOR(oldpi, nclass);
  MAKE_MATRIX(oldMu, nclass, p);
  MAKE_MATRIX(oldLTSigma, nclass, n_par);

  do{
/* Modified by Wei-Chen Chen on 2009/03/08.
    randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma);
    mb_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma);
*/
    randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma);

    iter = maxshortiter - totiter;
    if(fixed_iter > iter) fixed_iter = iter;
    iter = mod_shortemcluster(n, p, nclass, oldpi, X, oldMu, oldLTSigma,
                              fixed_iter, &llhval, conv_iter, conv_eps);
    if (llhval >= oldllh) {
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
} /* End of mod_shortems(). */
    

/* This function is equal to em_EM() with shorteps = Inf if fixed_iter = 1.
   This is a version for C.
*/
int rand_EM(double **x, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma, double *llhdval, int *nc, int shortiter, int fixediter,
    int *conv_iter, double *conv_eps){
  if(nclass == 1){
    meandispersion_MLE(x, n, p, Mu[0], LTSigma[0]);
    *llhdval = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
               0.5 * n * p * log(2 * PI);
  } else {
    mod_shortems(n, p, nclass, pi, x, Mu, LTSigma, shortiter, fixediter,
                 conv_iter, conv_eps);
    emcluster(n, p, nclass, pi, x, Mu, LTSigma, 1000, 0.0001, llhdval,
              conv_iter, conv_eps);
  } 

  return 0;
} /* End of rand_EM(). */

