/* This file contains several functions modified from Dr. Maitra's original
   codes in "init_rand_EM.c" and "M_init_rand_EM.c" to perform model-based
   initializers with the rand.EM procedure.
   Created: Wei-Chen Chen on 2009/02/07.
*/

#include "mb_tool.h"

/* Modified shortems() for model-based initializer. */
void mod_shortems_mb(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter){
  int i, iter, totiter = 0, n_par = p * (p + 1) / 2;
  double *oldpi, **oldMu, **oldLTSigma, oldllh = -Inf, llhval;

  MAKE_VECTOR(oldpi, nclass);
  MAKE_MATRIX(oldMu, nclass, p);
  MAKE_MATRIX(oldLTSigma, nclass, n_par);

  do{
    mb_init(X, n, p, nclass, oldpi, oldMu, oldLTSigma);

    iter = maxshortiter - totiter;
    if(fixed_iter > iter) fixed_iter = iter;
    iter = mod_shortemcluster(n, p, nclass, oldpi, X, oldMu, oldLTSigma,
                              fixed_iter, &llhval);
    
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
} /* End of mod_shortems_mb(). */
    


/* This is the model-based rand.EM.
   This function is equal to mb_em_EM() with shorteps = Inf if fixed_iter = 1.
*/
int mb_rand_EM(double **x, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma, double *llhdval, int shortiter, int fixed_iter){
  if(nclass == 1){
    meandispersion_MLE(x, n, p, Mu[0], LTSigma[0]);
    *llhdval = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
               0.5 * n * p * log(2 * PI);
  } else {
    mod_shortems_mb(n, p, nclass, pi, x, Mu, LTSigma, shortiter, fixed_iter);
    emcluster(n, p, nclass, pi, x, Mu, LTSigma, 1000, 0.0001, llhdval);
  } 

  return 0;
} /* End of mb_rand_EM(). */

