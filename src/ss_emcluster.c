/* This file is modified from "emcluster.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"

/* Modified by Wei-Chen Chen on 2009/02/02. */
void ss_estep(int n, int p, int k, double **X, double **Gamma, double *pi,
    double **Mu, double **LTSigma, int *lab){
  /* This is the E-Step: assigns the responsibilities of each of the k groups 
     to each of the n observations. The input are the n p-variate observations
     in the matrix of observations X to one of k multi-variate Gaussian 
     densities. Each density has prior class probabilities given by pi, and 
     has mean given by the kxp-matrix Mu and the kxp(p+1)/2 matrix of 
     corresponding lower triangular dispersion matrix vectors */

  /* Gamma[l][i]
     = pi_i * L_i / sum_j(pi_j * L_j)
     = 1 / sum_j(pi_j * L_j / (pi_i * L_i))
     = 1 / sum_j(exp(log(pi_j) + log(L_j) - log(pi_i) - log(L_i))) */
  int i, j, l, k1 = k - 1;
  double sum, sum_Gamma, log_pi[k], tmp_Gamma[k], *pt_Gamma;

  for(i = 0; i < k; i++){
    log_pi[i] = log(pi[i]);
  }

  for(l = 0; l < n; l++){
    if(lab[l] == -1){
      for(i = 0; i < k; i++){
        tmp_Gamma[i] = log_pi[i] + dlmvnorm(X[l], p, Mu[i], LTSigma[i]);
      }

      sum_Gamma = 0.0;
      pt_Gamma = Gamma[l];
      for (i = 0; i < k1; i++){
        sum = 0.0;
        for(j = 0; j < k; j++){
          sum = sum + exp(tmp_Gamma[j] - tmp_Gamma[i]);
        }
        sum = 1.0 / sum;
        sum_Gamma = sum_Gamma + sum; 
        /* Gamma[l][i] = sum; */
        *pt_Gamma = sum;
        pt_Gamma++;
      }
      /* Gamma[l][k1] = 1.0 - sum_Gamma; */
      *pt_Gamma = 1.0 - sum_Gamma;
    } else{
      for(i = 0; i < k; i++){
        Gamma[l][i] = (i != lab[l]) ? 0.0 : 1.0;
      }
    }
  }
  return;
} /* End of ss_estep(). */

void ss_emcluster_org(int n, int p, int k, double *pi, double **X, double **Mu, 
    double **LTSigma, int maxiter, double eps, double *llhdval, int *lab){
  int iter;
  double **gamm, llhd, oldllhd;

  MAKE_MATRIX(gamm, n, k);
  llhd = lnlikelihood(n, p, k, pi, X, Mu, LTSigma);
  iter = 0;
  do{
    oldllhd = llhd;
    ss_estep(n, p, k, X, gamm, pi, Mu, LTSigma, lab);
    mstep(X, n, p, k, pi, Mu, LTSigma, gamm);
    llhd = lnlikelihood(n, p, k, pi, X, Mu, LTSigma);
    iter++;
  } while((((oldllhd - llhd) / oldllhd) > eps) && (iter < maxiter));
  *llhdval = llhd;
  FREE_MATRIX(gamm);
  return;
} /* End of ss_emcluster_org(). */

