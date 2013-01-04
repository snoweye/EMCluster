/* This file is created by Wei-Chen Chen on 2009/06/15.
   Change EM to ME steps and ignore compute logL twice.
*/

#include<stdlib.h>
#include "array.h"
#include<math.h>
#include "mat_vec.h"

double lnlikelihood(int n, int p, int k, double *pi, double **X, double **Mu,
    double **LTSigma);
void mstep(double **X, int n, int p, int k, double *pi, double **Mu,
    double **LTSigma, double **Gamma);
void estep_gamma(int n, int p, int k, double **X, double **Gamma,
    double **Mu, double **LTSigma);
void norm_gamma(int n, int k, double **Gamma, double *pi);
double lnlikelihood_gamma(int n, int k, double **Gamma, double *pi);


int shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int maxiter, double eps, double *llhdval){
  int iter, i, n_par =  p * (p + 1) / 2;
  double *backup_pi, **backup_Mu, **backup_LTSigma;
  double **gamm, llhd, oldllhd, llh0;

  MAKE_VECTOR(backup_pi, k);
  MAKE_MATRIX(backup_Mu, k, p);
  MAKE_MATRIX(backup_LTSigma, k, n_par);
  MAKE_MATRIX(gamm, n, k);

  estep_gamma(n, p, k, X, gamm, Mu, LTSigma);
  llhd = lnlikelihood_gamma(n, k, gamm, pi);
  llh0 = llhd;
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
      break;
    }

    iter++;
  } while((fabs((oldllhd - llhd) / (llh0 - llhd)) > eps) && (iter < maxiter));
  *llhdval = llhd;

  FREE_VECTOR(backup_pi);
  FREE_MATRIX(backup_Mu);
  FREE_MATRIX(backup_LTSigma);
  FREE_MATRIX(gamm);
  return iter;
}

