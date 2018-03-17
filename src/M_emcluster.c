/* This file is created by Wei-Chen Chen on 2009/06/15.
   Change EM to ME steps and ignore compute logL twice.
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"
#include "mat_vec.h"


double dlmvnorm(double *x, int p, double *mu, double *LTsigma);
double lnlikelihood(int n, int p, int k, double *pi, double **X, double **Mu,
    double **LTSigma);
void mstep(double **X, int n, int p, int k, double *pi, double **Mu,
    double **LTSigma, double **Gamma);


/* Created by Wei-Chen Chen on 2009/06/15. */
double lnlikelihood_gamma(int n, int k, double **Gamma, double *pi){
  int i, l;
  double llhdval = 0, tmp;
  
  for(l = 0; l < n; l++){
    tmp = 0;
    for(i = 0; i < k; i++){
      tmp = tmp + pi[i] * exp(Gamma[l][i]);
    }
    llhdval = llhdval + log(tmp);
  }
  return llhdval;
} /* End of lnlikelihood_gamma(). */


/* Modified by Wei-Chen Chen on 2009/06/15. */
void estep_gamma(int n, int p, int k, double **X, double **Gamma,
    double **Mu, double **LTSigma){
  int i, l;

  for(l = 0; l < n; l++){
    for(i = 0; i < k; i++){
      Gamma[l][i] = dlmvnorm(X[l], p, Mu[i], LTSigma[i]);
    }
  }
  return;
} /* End of estep_gamma(). */


/* Modified by Wei-Chen Chen on 2009/06/17. */
void estep_unnorm_dlmvn(int n, int p, int k, double **X, double **Gamma,
    double *pi, double **Mu, double **LTSigma){
  int i, l;
  double logpi[k];

  for(i = 0; i < k; i++){
    logpi[i] = log(pi[i]);
  }

  for(l = 0; l < n; l++){
    for(i = 0; i < k; i++){
      Gamma[l][i] = logpi[i] + dlmvnorm(X[l], p, Mu[i], LTSigma[i]);
    }
  }
  return;
} /* End of estep_unnorm_dlmvn(). */


/* Modified by Wei-Chen Chen on 2009/06/15. */
void norm_gamma(int n, int k, double **Gamma, double *pi){
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
    for(i = 0; i < k; i++){
      tmp_Gamma[i] = log_pi[i] + Gamma[l][i];
    }

    sum_Gamma = 0.0;
    pt_Gamma = Gamma[l];
    for(i = 0; i < k1; i++){
      sum = 0.0;
      for(j = 0; j < k; j++){
        sum = sum + exp(tmp_Gamma[j] - tmp_Gamma[i]);
      }
      sum = 1.0 / sum;
      sum_Gamma = sum_Gamma + sum; 
      *pt_Gamma = sum;
      pt_Gamma++;
    }
    *pt_Gamma = 1.0 - sum_Gamma;
  }
  return;
} /* End of norm_gamma(). */


/* Modified by Wei-Chen Chen on 2009/06/15. */
void emcluster(int n, int p, int k, double *pi, double **X, double **Mu, 
    double **LTSigma, int maxiter, double eps, double *llhdval,
    int *conv_iter, double *conv_eps){
  int iter, i, n_par =  p * (p + 1) / 2;
  double *backup_pi, **backup_Mu, **backup_LTSigma;
  double **gamm, llhd, oldllhd;

  MAKE_VECTOR(backup_pi, k);
  MAKE_MATRIX(backup_Mu, k, p);
  MAKE_MATRIX(backup_LTSigma, k, n_par);
  MAKE_MATRIX(gamm, n, k);

  estep_gamma(n, p, k, X, gamm, Mu, LTSigma);
  llhd = lnlikelihood_gamma(n, k, gamm, pi);
  iter = 0;
  do{
    oldllhd = llhd;
    norm_gamma(n, k, gamm, pi);

    for(i = 0; i < k; i++) backup_pi[i] = pi[i];
    cpy(Mu, k, p, backup_Mu);
    cpy(LTSigma, k, n_par, backup_LTSigma);

    mstep(X, n, p, k, pi, Mu, LTSigma, gamm);
    estep_gamma(n, p, k, X, gamm, Mu, LTSigma);
    llhd = lnlikelihood_gamma(n, k, gamm, pi);

    if(oldllhd > llhd){
      for(i = 0; i < k; i++) pi[i] = backup_pi[i];
      cpy(backup_Mu, k, p, Mu);
      cpy(backup_LTSigma, k, n_par, LTSigma);
      llhd = oldllhd;
      iter--;
      break;
    }

    iter++;
    *conv_eps = fabs((oldllhd - llhd) / oldllhd);
  } while((*conv_eps > eps) && (iter < maxiter));
  *llhdval = llhd;
  *conv_iter = iter;

  FREE_VECTOR(backup_pi);
  FREE_MATRIX(backup_Mu);
  FREE_MATRIX(backup_LTSigma);
  FREE_MATRIX(gamm);
  return;
}

