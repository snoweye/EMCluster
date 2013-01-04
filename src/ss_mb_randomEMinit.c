/* This file is modified from "mb_randomEMinit.c" for semi-supervised
   clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"

int ss_mb_assign_closest(double **X, int n, int p, int nclass, double **Mu,
    double **LTSigma, int* clas, int *lab){
  int i, j, tmp_class;
  double tmp_pi, pi[nclass], **gamm, tmp_gamm;

  tmp_pi = 1.0 / (double) nclass;
  for(i = 0; i < nclass; i++) pi[i] = tmp_pi;
  
  /* Run one estep to update clas. */
  MAKE_MATRIX(gamm, n, nclass);
  ss_estep(n, p, nclass, X, gamm, pi, Mu, LTSigma, lab);
  for(i = 0; i < n; i++){
    if(lab[i] == -1){
      tmp_gamm = 0.0;
      tmp_class = 0;
      for(j = 0; j < nclass; j++){
        if(tmp_gamm < gamm[i][j]){
          tmp_gamm = gamm[i][j];
          tmp_class = j;
        }
      }
      clas[i] = tmp_class;
    } else{
      clas[i] = lab[i];
    }
  }

  FREE_MATRIX(gamm);
  return 0;  
} /* End of ss_mb_assign_closest(). */


/* This function is called by ss_shortems().
   Mu[0, ..., labK-1] should be assigned before calling this function.
*/
void ss_mb_randomEMinit(double **x, int n, int p, int nclass, double *pi,
    double **Mu, double **LTSigma,
    int *lab, int labK, int nonlab_total, int *lab_index){
  int *ordr, i, j, *clas, *nc;
  int n_par = p * (p + 1) / 2, tmp_j;
  int new_nclass = nclass - labK;
  double mu[p], ltsigma[n_par], tmp_sigma = 0.0, oldLTSigma[nclass][n_par];
  double labMu[labK][p];

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
  for(i = 0; i < labK; i++){
    for(j = 0; j < p; j++) labMu[i][j] = Mu[i][j];
  }
  
  /* Initial centers for all other unknown clusters. */
  MAKE_VECTOR(ordr, new_nclass);
  MAKE_VECTOR(clas, n);
  MAKE_VECTOR(nc, nclass);
  do{
    for(i = 0; i < labK; i++){
      for(j = 0; j < p; j++) Mu[i][j] = labMu[i][j];
    }
    i = srswor(nonlab_total, new_nclass, ordr);
    for(i = labK; i < nclass; i++){
      for(j = 0; j < p; j++) Mu[i][j] = x[lab_index[ordr[i - labK]]][j];
    }
    for(i = 0; i < nclass; i++){
      for(j = 0; j < n_par; j++) LTSigma[i][j] = oldLTSigma[i][j];
    }
    ss_mb_assign_closest(x, n, p, nclass, Mu, LTSigma, clas, lab);
    j = initials(x, n, p, nclass, nc, Mu, LTSigma, clas);
  } while(j == 0);
  for(i = 0; i < nclass; i++) pi[i] = 1. * nc[i] / n;
  FREE_VECTOR(nc);
  FREE_VECTOR(clas);
  FREE_VECTOR(ordr);
} /* End of ss_mb_randomEMinit(). */

