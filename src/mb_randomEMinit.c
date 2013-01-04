/* This file contains several functions to perform model-based initializers.
   Created: Wei-Chen Chen on 2009/03/04.
*/

#include "mb_tool.h"

int mb_assign_closest(double **X, int n, int p, int nclass, double **Mu,
    double **LTSigma, int* clas){
  int i, j, tmp_class;
  double tmp_pi, pi[nclass], **gamm, tmp_gamm;

  tmp_pi = 1.0 / (double) nclass;
  for(i = 0; i < nclass; i++) pi[i] = tmp_pi;
  
  /* Run one estep to update clas. */
  MAKE_MATRIX(gamm, n, nclass);
  estep(n, p, nclass, X, gamm, pi, Mu, LTSigma);
  for(i = 0; i < n; i++){
    tmp_gamm = 0.0;
    tmp_class = 0;
    for(j = 0; j < nclass; j++){
      if(tmp_gamm < gamm[i][j]){
        tmp_gamm = gamm[i][j];
        tmp_class = j;
      }
    }
    clas[i] = tmp_class;
  }

  FREE_MATRIX(gamm);
  return 0;  
} /* End of mb_assign_closest(). */

void mb_randomEMinit(double **x, int n, int p, int nclass, double *pi,
    double **Mu, double **LTSigma){
  int *ordr, i, j, *clas, *nc;
  int n_par = p * (p + 1) / 2, tmp_j;
  double mu[p], ltsigma[n_par], tmp_sigma = 0.0, oldLTSigma[nclass][n_par];

  /* Initial a common sphere structure for all clusters. */
  meandispersion_MLE(x, n, p, mu, ltsigma);
  tmp_sigma = ltsigma[0];
  for(j = 1; j < p; j++){
    if(tmp_sigma < ltsigma[j]) tmp_sigma = ltsigma[j];
  }
  for(i = 0; i < nclass; i++){
    for(j = 0; j < n_par; j++) oldLTSigma[i][j] = 0.0;
    tmp_j = 0;
    for(j = 0; j < p; j++){
      oldLTSigma[i][tmp_j] = tmp_sigma;
      tmp_j = tmp_j + p - j;
    }
  }
  
  /* Initial centers for all clusters. */
  MAKE_VECTOR(ordr, nclass);
  MAKE_VECTOR(clas, n);
  MAKE_VECTOR(nc, nclass);
  do{
    i = srswor(n, nclass, ordr);
    for(i = 0; i < nclass; i++){
      for(j = 0; j < p; j++) Mu[i][j] = x[ordr[i]][j];
    }
    for(i = 0; i < nclass; i++){
      for(j = 0; j < n_par; j++) LTSigma[i][j] = oldLTSigma[i][j];
    }
    mb_assign_closest(x, n, p, nclass, Mu, LTSigma, clas);
    j = initials(x, n, p, nclass, nc, Mu, LTSigma, clas);
  } while(j == 0);
  for(i = 0; i < nclass; i++) pi[i] = 1. * nc[i] / n;
  FREE_VECTOR(nc);
  FREE_VECTOR(clas);
  FREE_VECTOR(ordr);
} /* End of mb_randomEMinit(). */


