#include<stdlib.h>
#include<stdio.h>
#include "array.h"
#include<math.h>
//WCC #define PI 3.141593
#define Inf 1e+140

#include <R.h>
#include <Rmath.h>

double dlmvnorm(double *x, int p, double *mu, double *LTsigma);
double dlmvnorm_singular(double *x, int p, double *mu, double *LTsigma);

int classify(double *X,int p,int k,double *pi, double **Mu, double **LTSigma)
{
  int j,l,class=0;
  double *mu,*ltsigma,temp,dum,dum1;

  MAKE_VECTOR(mu,p);
  MAKE_VECTOR(ltsigma,p*(p+1)/2);
  temp=-Inf;
  for (l=0;l<k;l++) {
    for (j=0;j<p;j++) {
      mu[j]=Mu[l][j];
    }
    for (j=0;j<(p*(p+1)/2);j++) {
      ltsigma[j]=LTSigma[l][j];
    }
    dum1=dlmvnorm(X,p,mu,ltsigma);
    /*dum1=dlmvnorm_singular(X,p,mu,ltsigma);*/
    dum=log(pi[l])+dum1;
    if (dum>temp) {
      temp=dum;
      class=l;
    }
  }
  FREE_VECTOR(mu);
  FREE_VECTOR(ltsigma);
  return(class);
}

