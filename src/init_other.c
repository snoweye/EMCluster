#include<stdlib.h>
#include "array.h"
#include<math.h>
#include "mat_vec.h"
// #define PI 3.141593
#define Inf 1e+140

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

int srswor(int n, int k, int *ranordr);
int assign_closest(double *X,int p,int nclass,double **Mu);
int initials(double **x,int n,int p,int nclass,int *nc,
	     double **Mu,double **LTSigma,int *class);
int emcluster(int n,int p,int k,double *pi,double **X,double **Mu, 
	      double **LTSigma,int maxiter,double eps,double *llhdval,
	      int *conv_iter,double *conv_eps);
double determinant(double *LTSigma,int n);
void assign(int n, int p,int k,double **X,double *pi,double **Mu,
            double **LTSigma,int *class,int *nc);
/* Modified: Wei-Chen Chen on 2008/12/05.
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);
*/
void meandispersion_MLE(double **x, int n, int p, double *mu, double *ltsigma);
double aic(double llhd,int nobs, int ndim, int nclus,int aicbic);
double lnlikelihood(int n,int p,int k,double *pi,double **X,double **Mu,
                    double **LTSigma);
void estep(int n,int p,int k,double **X,double **Gamma,double *pi,double **Mu, 
           double **LTSigma);
void mstep(double **X,int n,int p,int k,double *pi,double **Mu,
           double **LTSigma,double **Gamma);

/* Add by Wei-Chen Chen on 2009/03/08. */
int mb_randomEMinit(double **x, int n, int p, int nclass, double *pi,
    double **Mu, double **LTSigma);
int shortemcluster(int n,int p,int k,double *pi,double **X,double **Mu,  
		   double **LTSigma,int maxiter,double eps,double *llhdval,
		   int *conv_iter,double *conv_eps);

int randomEMinit(double **x, int n, int p, int nclass, double *pi,
		 double **Mu, double **LTSigma)
{
  int *ordr,i,j,*clas,*nc;
  /* This is a bug. Modified by Wei-Chen Chen on 2009/03/13.
    MAKE_VECTOR(ordr,n);
  */
  MAKE_VECTOR(ordr, nclass);
  MAKE_VECTOR(clas, n);
  MAKE_VECTOR(nc, nclass);
  do {
    /* This is a bug. Modified by Wei-Chen Chen on 2009/03/13.
      i=srswor(n, n, ordr);
    */
    i=srswor(n, nclass, ordr);
    for(i=0;i<nclass;i++){
      for(j=0;j<p;j++) Mu[i][j]=x[ordr[i]][j];
    }
    for(i=0;i<n;i++) clas[i]=assign_closest(x[i],p,nclass,Mu);
    j=initials(x,n,p,nclass,nc,Mu,LTSigma,clas);
  } while (j==0);
  for(i=0;i<nclass;i++) pi[i]=1.*nc[i]/n;
  FREE_VECTOR(nc);
  FREE_VECTOR(clas);
  FREE_VECTOR(ordr);
  return 0;
}


int shortemcluster_org(int n,int p,int k,double *pi,double **X,double **Mu,  
		   double **LTSigma,int maxiter,double eps,double *llhdval)
{
  int iter;
  double **gamm,llhd,oldllhd,llh0;
  /*same as emcluster, only difference being in how convergence is handled,
   done as per Biernacki et al*/
  MAKE_MATRIX(gamm,n,k);
  llhd=lnlikelihood(n,p,k,pi,X,Mu,LTSigma);
  llh0=llhd;
  iter=0;
  do  {
    oldllhd=llhd;
    estep(n,p,k,X,gamm,pi,Mu,LTSigma);
    mstep(X,n,p,k,pi,Mu,LTSigma,gamm);
    llhd=lnlikelihood(n,p,k,pi,X,Mu,LTSigma);
    iter++;
  } while (((oldllhd-llhd)/(llh0-llhd) > eps) && (iter<maxiter));
/* Marked by Wei-Chen Chen on 2009/01/21.
   This value, (oldllhd-llhd)/(llhd-llh0), is always negative.

  } while ((((oldllhd-llhd)/(llhd-llh0)) > eps) && (iter<maxiter));
*/
  (*llhdval)=llhd;
  FREE_MATRIX(gamm);
  return iter;
}


int shortems(int n,int p,int nclass,double *pi,double **X,double **Mu,  
             double **LTSigma,int maxshortiter,double shorteps,
             int *conv_iter, double *conv_eps)
{
  /*initializing as per Beiernacki, Celeaux, Govaert~(2003) */

  int i,iter,totiter=0;
  double *oldpi,**oldMu,**oldLTSigma,oldllh=-Inf,llhval;
  MAKE_VECTOR(oldpi,nclass);
  MAKE_MATRIX(oldMu,nclass,p);
  MAKE_MATRIX(oldLTSigma,nclass,p*(p+1)/2);
  do {
/* Modified by Wei-Chen Chen on 2009/03/08.
    i=randomEMinit(X,n,p,nclass,oldpi,oldMu,oldLTSigma);
    i = mb_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma);
*/
    i = randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma);

    iter=maxshortiter-totiter;
    iter=shortemcluster(n,p,nclass,oldpi,X,oldMu,oldLTSigma,iter,shorteps,
			&llhval,conv_iter,conv_eps);
    if (llhval >= oldllh) {
      int i;
      oldllh=llhval;
      cpy(oldMu,nclass,p,Mu);
      cpy(oldLTSigma,nclass,p*(p+1)/2,LTSigma);
      for(i=0;i<nclass;i++) pi[i]=oldpi[i];
    }
    totiter+=iter;
  }  while (totiter < maxshortiter);
  FREE_MATRIX(oldMu);
  FREE_MATRIX(oldLTSigma);
  FREE_VECTOR(oldpi);
  return totiter; 
}
    

/*
  This routine clusters a dataset into an optimal number of groups using the 
  E-M algorithm and a choice of AIC or BIC and with start points provided by 
  the shortems function.  
  The function returns    the classification ids as well as the parameter 
  estimates. The parameter estimates are created in the routine and should not
  be allocated before the function call. This also means that the parameter 
  estimates have to be cleaned after the function call. (NOTE: This is an 
  important fact to remember. */
int em_EM(double **x,int n,int p,int nclass,double *pi,double **Mu,
	  double **LTSigma,double *llhdval,int *nc,int shortiter,
	  double shorteps,int *conv_iter,double *conv_eps)
{
  int flag=0;
  double like;

  if (nclass==1) {
/* These formulae are not correct with two errors,
   1. n-1 should be replace by n in meandispersion() for LTSigma, and
   2. determinant() will change the values of LTSigma[0], see "initials.c".
   Modified: Wei-Chen Chen on 2008/12/05.

    meandispersion(x,n,p,Mu[0],LTSigma[0]);
*/
    meandispersion_MLE(x,n,p,Mu[0],LTSigma[0]);
    like=-0.5*n*p-0.5*n*log(determinant(LTSigma[0],p))-0.5*n*p*log(2*M_PI);
    (*llhdval)=like;
  }
  else {
    shortems(n,p,nclass,pi,x,Mu,LTSigma,shortiter,shorteps,conv_iter,conv_eps);
    emcluster(n,p,nclass,pi,x,Mu,LTSigma,1000,0.0001,&like,conv_iter,conv_eps);
    (*llhdval)=like;
  } 
  return flag;
}


void runemEMcluster(double **x,int n,int p,int kmin,int kmax,int *kopt,
		    double *pi,double **Mu,double **LTSigma,int *nc,
		    int *class,int aicind,double alpha,int shortiter,
		    double shorteps,int *conv_iter,double *conv_eps)
{
  int iclass,j,*currnc;
  double **currMu,**currLTSigma,*currpi,aicval,currval=Inf,llhd;
  for(iclass=kmin;iclass<=kmax;iclass++) {
    MAKE_MATRIX(currMu,iclass,p);
    MAKE_MATRIX(currLTSigma,iclass,p*(p+1)/2);
    MAKE_VECTOR(currpi,iclass);
    MAKE_VECTOR(currnc,iclass);
    if(!em_EM(x,n,p,iclass,currpi,currMu,currLTSigma,&llhd,currnc,shortiter,shorteps,conv_iter,conv_eps)) {
      aicval=aic(llhd,n,p,iclass,aicind);
      if (aicind==1) 
        Rprintf("nclass= %d like = %f AIC = %f \n",iclass,llhd,aicval);
      else 
        Rprintf("nclass = %d like = %f BIC = %f \n",iclass,llhd,aicval);
      if(aicval<currval) {
        currval=aicval;
        (*kopt)=iclass;
        j=cpy(currMu,iclass,p,Mu);
        j=cpy(currLTSigma,iclass,p*(p+1)/2,LTSigma);
        for (j=0;j<iclass;j++) {
          nc[j]=currnc[j];
          pi[j]=currpi[j];
        }
      }
    }
    FREE_MATRIX(currMu);
    FREE_MATRIX(currLTSigma);
    FREE_VECTOR(currpi);
    FREE_VECTOR(currnc);
  }
  if (*kopt==1) for(j=0;j<n;j++) class[j]=0;
  else assign(n,p,(*kopt),x,pi,Mu,LTSigma,class,nc);
  return;
}

