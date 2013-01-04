/*
  This routine clusters a dataset into groups using the E-M algorithm with 
  start points provided by the starts_via_svd routine. The function returns 
  the classification ids as well as the parameter estimates.
*/

#include<stdlib.h>
#include<stdio.h>
#include "array.h"
#include<math.h>
// #define PI 3.141593
#define Inf 1e+140

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double determinant(double *LTSigma,int n);
void emcluster(int n,int p,int k,double *pi,double **X,double **Mu, 
	       double **LTSigma,int maxiter,double eps,double *llhdval);
int starts_via_svd(int n,int m,double **Mu,double **x,int nclus,int *ningrp,
                    double *pi,int *grpids,double **LTSigma,double alpha,
                    int llhdnotW);
void assign(int n, int p,int k,double **X,double *pi,double **Mu,
	    double **LTSigma,int *class,int *nc);
/* Modified: Wei-Chen Chen on 2008/12/05.
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);
*/
void meandispersion_MLE(double **x, int n, int p, double *mu, double *ltsigma);
int emgroup(double **x,int n,int p,int nclass,double *pi,double **Mu,
	     double **LTSigma,double *llhdval,int *nc,int *class)
{
  int j,flag=0;
  double like;
  
  if (nclass==1) {
    nc[nclass-1]=n;
    pi[nclass-1]=1.0;
    for (j=0;j<n;j++) class[j]=0;
/* These formulae are not correct with two errors,
   1. n-1 should be replace by n in meandispersion() for LTSigma, and
   2. determinant() will change the values of LTSigma[0], see "initials.c".
   Modified: Wei-Chen Chen on 2008/12/05.
    meandispersion(x,n,p,Mu[0],LTSigma[0]);
*/
    meandispersion_MLE(x,n,p,Mu[0],LTSigma[0]);
    like=-0.5*n*p-0.5*n*log(determinant(LTSigma[0],p))-0.5*n*p*log(2*PI);
  }
  else {
    if(!starts_via_svd(n,p,Mu,x,nclass,nc,pi,class,LTSigma,0.99,1))      {
      for(j=0;j<nclass;j++) pi[j]=nc[j]/(double)n;
      emcluster(n,p,nclass,pi,x,Mu,LTSigma,1000,0.0001,&like);
      assign(n,p,nclass,x,pi,Mu,LTSigma,class,nc);
    }
    else flag=1;
  }
  (*llhdval)=like;
  Rprintf("like =  %f\n",like);
  return flag;
} 

