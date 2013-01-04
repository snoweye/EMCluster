#include<stdlib.h>
#include<stdio.h>
#include "array.h"
#include "order.h"
#include<math.h>

//WCC #define PI 3.141593
#define Inf 1e+140

#include <R.h>
#include <Rmath.h>

int pposymatinv(int N,double *A, char UPLO, double *determinant);
void hclass(int n, int *ia, int *ib, int lev, int *iclass);
int classify(double *X,int p,int k,double *pi, double **Mu, double **LTSigma);
void unique(int n,double *x,int *m,double *y);

double determinant(double *LTSigma,int n)
{
  double dum;

/* Bugs: pposymatinv() and posymatinv() in "inverse.c" will modify LTSigma.
   Make a copy before call pposymatinv().
   Modified: Wei-Chen Chen on 2008/12/03.

   pposymatinv(n,LTSigma,'U',&dum);
*/ 

  double *a;
  int i;
  MAKE_VECTOR(a,n*(n+1)/2);
  for (i=0;i<n*(n+1)/2;i++) a[i]=LTSigma[i];
  pposymatinv(n,a,'U',&dum);
  FREE_VECTOR(a);

  return dum;
}

void meandispersion(double **x, int n, int p, double *mu, double *ltsigma)
{
  /* This routine calculates the mean and dispersion of a homogeneous sample */
  int i,j,l;

  for (i=0;i<(p*(p+1)/2);i++) ltsigma[i]=0.;
  for (i=0;i<p;i++) mu[i]=0.;
  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) mu[j]+=x[i][j];
  }
  for (j=0;j<p;j++) mu[j]/=n;
  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) {
      for (l=0;l<=j;l++) ltsigma[j*(j+1)/2+l]+=(x[i][j]-mu[j])*(x[i][l]-mu[l]);
    }
  }
  if(n>1) {
    for (j=0;j<(p*(p+1)/2);j++) ltsigma[j]/=n-1;
  }
  return;
}


int initials(double **x,int n,int p,int nclass,int *nc,
	      double **Mu,double **LTSigma,int *class)
{
  double **y;
  int i,j,k,l,m=1;
  for (i=0;i<nclass;i++) {
       nc[i]=0;
       for(l=0;l<n;l++) if (class[l]==i) nc[i]++;
  }
  for(i=0;i<nclass;i++) {
    if(nc[i]>p) m*=1;
    else m*=0;
    MAKE_MATRIX(y,nc[i],p);
    k=0;
    for(l=0;l<n;l++) {
      if (class[l]==i) {
	for (j=0;j<p;j++)   y[k][j]=x[l][j];
	k++;
      }
    }
    meandispersion(y,nc[i],p,Mu[i],LTSigma[i]);
    FREE_MATRIX(y);
  }
  return m;
}

void assign(int n, int p,int k,double **X,double *pi,double **Mu,
	   double **LTSigma,int *class,int *nc) 
{
  int i,j;
  double *x;

  MAKE_VECTOR(x,p);
  for (i=0;i<k;i++)  nc[i]=0;
  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) x[j]=X[i][j];
    class[i]=classify(x,p,k,pi,Mu,LTSigma);
    nc[class[i]]++;
  }
  FREE_VECTOR(x);
  return;
}

double aic(double llhd,int nobs, int ndim, int nclus,int aicbic)
{
  int noparameters;
  noparameters=ndim*nclus+nclus*ndim*(ndim+1)/2+(nclus-1);
  if (aicbic==1) {
    return(-2*llhd+2.0*noparameters); /*This is  the AIC*/
  }
  else {
  return(-2*llhd+noparameters*log((double)nobs)); /*This is the BIC*/
  }
}

void starters(double **x,int n,int p,int nclass,double **Mu,double *pi,
	      double **LTSigma,int *ia,int *ib,double *crit,int *iflag)
/* This routine takes the output from hc and cuts the tree into the first 
   nclass groups giving nclass nonsingular estimates of dispersions. This is 
   achieved by recursively looking at the smallest k>=nclass which has exactly
   nclass many groups with number of observations > p. If there such a k,
   initial estimates of Mu, LTSigma, and pi, based only on estimates from 
   clusters with number of observations > p are passed on for use as initial
   estimates in the emcluster routine. If there is no such k, (k=n) and iflag=1
   and the other estimates are meaningless and should be ignored.*/ 
{
  int m,i,j,k=nclass,*nc,*class,numpd,tag=0,kk;
  
  MAKE_VECTOR(nc,n);
  MAKE_VECTOR(class,n);
  (*iflag)=0;
  while ((k<n) && (tag==0)) {
    hclass(n,ia,ib,k,class);
    for (i=0;i<k;i++) {
      nc[i]=0;
    }
    for (i=0;i<n;i++) {
      nc[class[i]]++;
    }
    numpd=0;
    for (i=0;i<k;i++) {
      if (nc[i]>p) {
	numpd++;
      }
    }
    if (numpd>=nclass) {
      double **y;
      MAKE_MATRIX(y,n,p);
      tag=1;
      if (k==nclass) {
	m=n;
	for (i=0;i<m;i++) {
	  for (j=0;j<p;j++) {
	    y[i][j]=x[i][j];
	  }
	}
      }
      else {
	double *mc,*dum1;
	int *idx,*itag,*newclass;
	MAKE_VECTOR(mc,k); /* we need to get the nclass largest clusters*/
	MAKE_VECTOR(itag,k);
	MAKE_VECTOR(newclass,n);
	for (i=0;i<k;i++) {
	  itag[i]=0;
	}
	MAKE_VECTOR(mc,k);
	for (i=0;i<k;i++) {
	  mc[i]=-1.0*nc[i];
	}
	idx=(int *)orderDouble(mc,k);
	FREE_VECTOR(mc);
	for (i=0;i<nclass;i++) {
	  itag[idx[i]]=1;
	}
	
	m=0; 
	for (i=0;i<n;i++) {
	  if (itag[class[i]]==1) {
	    for (j=0;j<p;j++) {
	      y[m][j]=x[i][j];
	      newclass[m]=class[i];
	  }
	    m++;
	  }
	}
	FREE_VECTOR(itag);
	/* now we need to re-assign the class indicators by getting the unique 
	   sorted values and then sort them and assign them again*/
	MAKE_VECTOR(mc,m);
	MAKE_VECTOR(dum1,m);
	MAKE_VECTOR(itag,m);
	for (i=0;i<m;i++) {
	  mc[i]=newclass[i];
	}
	unique(m,mc,&kk,dum1);
	for (i=0;i<m;i++) {
	  itag[i]=0;
	}
	for (j=0;j<kk;j++) {
	  for (i=0;i<m;i++) {
	    if ((((int)dum1[j])==newclass[i]) && (itag[i]==0)) {
	      class[i]=j;
	      itag[i]=1;
	    }
	  }
	}
	FREE_VECTOR(dum1);
	FREE_VECTOR(mc);
	FREE_VECTOR(itag);
	FREE_VECTOR(idx);
	FREE_VECTOR(newclass);
      }
      
      for (i=0;i<nclass;i++) {
	nc[i]=0;
      }
      for (i=0;i<m;i++) {
	nc[class[i]]++;
      }
      initials(y,m,p,nclass,nc,Mu,LTSigma,class);
      for (i=0;i<nclass;i++) {
	pi[i]=nc[i]/(double)m;
      }
      FREE_MATRIX(y);
    }
    else {
      k++;
    }
  }
  if (k==n) {
    (*iflag)=1;
  }
  FREE_VECTOR(nc);
  FREE_VECTOR(class);
  return;
}
  
