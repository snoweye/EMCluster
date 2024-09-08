/* This file is modified from "rand_EM.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"

int ss_mod_shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int fixed_iter, double *llhdval,
    int *conv_iter, double *conv_eps, int *lab){
  int iter;
  double **gamm;

  MAKE_MATRIX(gamm,n,k);
  iter=0;
  do{
    ss_estep(n, p, k, X, gamm, pi, Mu, LTSigma, lab);
    mstep(X, n, p, k, pi, Mu, LTSigma, gamm);
    iter++;
  } while(iter < fixed_iter);
  *conv_iter = iter;
  *conv_eps = -1.0;

  *llhdval = lnlikelihood(n, p, k, pi, X, Mu, LTSigma);
  FREE_MATRIX(gamm);

  return iter;
} /* End of ss_mod_shortemcluster(). */


void ss_mod_shortems(int n, int p, int nclass, double *pi, double **X,
    double **Mu, double **LTSigma, int maxshortiter, int fixed_iter,
    int *conv_iter, double *conv_eps, int *lab, int labK){
  int i, j, iter, totiter = 0, n_par = p * (p + 1) / 2;
  int nonlab_total = 0, lab_index[n];
  double *oldpi, **oldMu, **oldLTSigma, oldllh = -Inf, llhval;
  double **labMu;

  MAKE_VECTOR(oldpi, nclass);
  MAKE_MATRIX(oldMu, nclass, p);
  MAKE_MATRIX(oldLTSigma, nclass, n_par);
  MAKE_MATRIX(labMu, labK, p);

  for(i = 0; i < n; i++){
    if(lab[i] == -1) lab_index[nonlab_total++] = i;
  }
  labInitMus(n, p, labK, X, lab, labMu);

  do{
    for(i = 0; i < labK; i++){
      for(j = 0; j < p; j++) oldMu[i][j] = labMu[i][j];
    }

    iter = maxshortiter - totiter;

/* Modified by Wei-Chen Chen on 2009/03/08.
    ss_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma,
                    lab, labK, nonlab_total, lab_index);
    ss_mb_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma,
                       lab, labK, nonlab_total, lab_index);
*/
    ss_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma,
                    lab, labK, nonlab_total, lab_index);

    if(fixed_iter > iter) fixed_iter = iter;
    iter = ss_mod_shortemcluster(n, p, nclass, oldpi, X, oldMu, oldLTSigma,
                                 fixed_iter, &llhval, conv_iter, conv_eps, lab);

    if (llhval >= oldllh) {
      oldllh = llhval;
      cpy(oldMu, nclass, p, Mu);
      cpy(oldLTSigma, nclass, n_par, LTSigma);
      for(i = 0; i < nclass; i++) pi[i] = oldpi[i];
    }

    totiter += iter;
  } while(totiter < maxshortiter);

  FREE_MATRIX(oldMu);
  FREE_MATRIX(oldLTSigma);
  FREE_VECTOR(oldpi);
  FREE_MATRIX(labMu);
} /* End of ss_mod_shortems(). */


/* This function is equal to em_EM() with shorteps = Inf if fixed_iter = 1.
   This is a version for C.
*/
int ss_rand_EM(double **x, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma, double *llhdval, int *nc, int shortiter, int fixediter,
    int *conv_iter, double *conv_eps, int *lab, int labK){
  if(nclass == 1){
    meandispersion_MLE(x, n, p, Mu[0], LTSigma[0]);
    *llhdval = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
               0.5 * n * p * log(2 * M_PI);
  } else {
    ss_mod_shortems(n, p, nclass, pi, x, Mu, LTSigma, shortiter, fixediter,
                    conv_iter, conv_eps, lab, labK);
    ss_emcluster(n, p, nclass, pi, x, Mu, LTSigma, 1000, 0.0001, llhdval,
                 conv_iter, conv_eps, lab);
  } 

  return 0;
} /* End of ss_rand_EM(). */

