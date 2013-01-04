#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"
//WCC #define PI 3.141593
#define Inf 1e+140

#include <R.h>
#include <Rmath.h>

double dlmvnorm(double *x, int p, double *mu, double *LTsigma);
double dlmvnorm_singular(double *x, int p, double *mu, double *LTsigma);
double lnlikelihood(int n,int p,int k,double *pi,double **X,double **Mu,
		    double **LTSigma);

void estep_org(int n,int p,int k,double **X,double **Gamma,double *pi,double **Mu, 
	   double **LTSigma)
{
  /* This is the E-Step: assigns the responsibilities of each of the k groups 
     to each of the n observations. The input are the n p-variate observations
     in the matrix of observations X to one of k multi-variate Gaussian 
     densities. Each density has prior class probabilities given by pi, and 
     has mean given by the kxp-matrix Mu and the kxp(p+1)/2 matrix of 
     corresponding lower triangular dispersion matrix vectors */

  int i,l;
  double temp,sum;

  for (l=0;l<n;l++) {
    sum=0.;
    for (i=0;i<k;i++) {
      temp=dlmvnorm(X[l],p,Mu[i],LTSigma[i]);
      /*      temp=dlmvnorm_singular(X[l],p,Mu[i],LTSigma[i]);*/
      Gamma[l][i]=pi[i]*exp(temp);
      sum+=pi[i]*exp(temp);
    }
    for (i=0;i<k;i++) {
      Gamma[l][i]/=sum;
    }
  }
  return;
}


/* Modified by Wei-Chen Chen on 2009/02/02. */
void estep(int n, int p, int k, double **X, double **Gamma, double *pi,
    double **Mu, double **LTSigma) {
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
/*      Gamma[l][i] = sum; */
      *pt_Gamma = sum;
      pt_Gamma++;
    }
/*    Gamma[l][k1] = 1.0 - sum_Gamma; */
    *pt_Gamma = 1.0 - sum_Gamma;
  }
  return;
} /* End of estep(). */

void mstep(double **X,int n,int p,int k,double *pi,double **Mu,
	   double **LTSigma,double **Gamma) 
{
  int i,j,l,ll;
  double *sum,sumpi;

  MAKE_VECTOR(sum,k);
  for (ll=0;ll<k;ll++) {
    sum[ll]=0.;
    for (j=0;j<p;j++) {
      Mu[ll][j]=0.;
    }
    for (j=0;j<(p*(p+1)/2);j++) {
      LTSigma[ll][j]=0.;
    }
    for (i=0;i<n;i++) {
      sum[ll]+=Gamma[i][ll];
      for (j=0;j<p;j++) {
	Mu[ll][j]+=Gamma[i][ll]*X[i][j];
      }
    }
    for (j=0;j<p;j++) {
      Mu[ll][j]/=sum[ll];
    }
    for (i=0;i<n;i++) {
      for (j=0;j<p;j++) {
	for (l=0;l<=j;l++) {
	  LTSigma[ll][j*(j+1)/2+l]+=Gamma[i][ll]*(X[i][j]-Mu[ll][j])*(X[i][l]-Mu[ll][l]);
	}
      }
    }
    for (j=0;j<p;j++) {
      for (l=0;l<=j;l++) {
        LTSigma[ll][j*(j+1)/2+l]/=sum[ll];
      }
    }
  }
  sumpi=0.;
  for (i=0;i<k;i++) {
    sumpi+=sum[i];
  }
  for (i=0;i<k;i++) {
    pi[i]=sum[i]/sumpi;
  }
  FREE_VECTOR(sum);
  return;
}

void emcluster_org(int n,int p,int k,double *pi,double **X,double **Mu, 
	     double **LTSigma,int maxiter,double eps,double *llhdval)
{
  int iter;
  double **gamm,llhd,oldllhd;
/*
*  for (i=0;i<k;i++) {
*    printf("%f ",pi[i]);
*  }
*  printf("\n");
*/
  MAKE_MATRIX(gamm,n,k);
  llhd=lnlikelihood(n,p,k,pi,X,Mu,LTSigma);
  iter=0;
  do  {
    oldllhd=llhd;
    estep(n,p,k,X,gamm,pi,Mu,LTSigma);
    mstep(X,n,p,k,pi,Mu,LTSigma,gamm);
    llhd=lnlikelihood(n,p,k,pi,X,Mu,LTSigma);
    iter++;
  } while ((((oldllhd-llhd)/oldllhd) > eps) && (iter<maxiter));
  (*llhdval)=llhd;
  FREE_MATRIX(gamm);
  return;
}

