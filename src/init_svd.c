#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_vec.h"
#include "array.h"
#include "order.h"
#include "quantile.h"
#define Inf 1e+140
#define SQ(x) ((x)*(x))

void hclassify(int n,int m, double **x,int hcrit,int nclass,int *class);

void prcomp(int n, int m,double **x,double **UtX,double *D);

void kmeans(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
	    int iter,double *wss,int *ifault);
int initials(double **x,int n,int p,int nclass,int *nc,
	     double **Mu,double **LTSigma,int *class);
double determinant(double *LTSigma,int n);
double  lnlikelihood(int n,int p,int k,double *pi,double **X,double **Mu,
		     double **LTSigma);
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);

int svdd(double **a, int m, int n, double *d, double **u, double **v);

double detw(int nclass,int p,double **LTSigma,int *nc)
{
  double *temp,detW;
  int i,j;
  MAKE_VECTOR(temp,p*(p+1)/2);
  for(j=0;j<(p*(p+1)/2);j++) temp[j]=0;
  for(i=0;i<nclass;i++) {
    for (j=0;j<(p*(p+1)/2);j++) temp[j]+=(nc[i]-1)*LTSigma[i][j];
  }
  detW=determinant(temp,p); 
  FREE_VECTOR(temp);
  return detW;
}

/* write the digits of n in its broken down "base" into buf */
void break_down(int n, int *base, int *buf, int buflen)
{
        int i;
        for (i=0; i<buflen; i++) {
	  buf[i] = n%base[i];
	  n /= base[i];
        }
}

int assign_closest(double *X,int p,int nclass,double **Mu)
{
  int j,l,class=0;
  double temp,dum1;

  temp=Inf;
  for (l=0;l<nclass;l++) {
    dum1=0.;
    for(j=0;j<p;j++) dum1+=SQ(X[j]-Mu[l][j]);
    if (dum1<temp) {
      temp=dum1;
      class=l;
    }
  }
  return class;
}

double **eliminulls(double **x,int n,int p,int *nclass,double **Mu,int kk)
{ 
  /*This routine eliminates all those centers which have representation less
    than or equal to kk and returns the remaining centers */

  int i,j,*nc,k=0,newtotcl=(*nclass);
  double **Centers;

  MAKE_VECTOR(nc,(*nclass));
  for(i=0;i<(*nclass);i++) nc[i]=0;
  for(i=0;i<n;i++) nc[assign_closest(x[i],p,(*nclass),Mu)]++;
  for(i=0;i<(*nclass);i++) if(nc[i]<=kk) newtotcl--;

  MAKE_MATRIX(Centers,newtotcl,p);
  for(i=0;i<(*nclass);i++){
    if (nc[i]>kk) {
      for(j=0;j<p;j++) Centers[k][j]=Mu[i][j];
      k++;
    }
  }
  (*nclass)=newtotcl;
  FREE_VECTOR(nc);
  return Centers;
}

int starts_in_svd_domain(int n,int m,double **y,int nclus,int dmin,double *D,
			  int *ningrp,int *grpids)
{    
  double **cent,*dum1,*probs,*qtl,*cctr,**dumy,**dumcctr,
    *wss1,**ccent,**mu,**ltsigma;
  int i,j,k,*counts,*ncl,totcl=1,sumcl=0,*buf,*dum2,*ning,kopt,
    *grclass,*cum1,maxiter=100000,ind=1;   

  MAKE_VECTOR(ncl,dmin);
  i=(int)ceil(pow((m-dmin+1)*nclus,1.0/dmin));
  for(j=0;j<dmin;j++) ncl[j]=i*round(D[j]/D[dmin-1]);
  for(j=0;j<dmin;j++) {
    totcl*=ncl[j];
    sumcl+=ncl[j];
  }
  
  MAKE_VECTOR(dum1,n);
  MAKE_VECTOR(cctr,sumcl);
  MAKE_MATRIX(dumy,n,1);
  MAKE_VECTOR(dum2,n);

  k=0;
  for(i=0;i<dmin;i++) {
    MAKE_VECTOR(ning,ncl[i]);
    MAKE_VECTOR(qtl,ncl[i]);
    MAKE_MATRIX(dumcctr,ncl[i],1);
    MAKE_VECTOR(wss1,ncl[i]);
    for(j=0;j<n;j++) dum1[j]=y[j][i];
    MAKE_VECTOR(probs,ncl[i]);
    for(j=0;j<ncl[i];j++) probs[j]=j/(ncl[i]-1.);
    j=quantile(n,dum1,probs,qtl,ncl[i]);
    for(j=0;j<n;j++) dumy[j][0]=y[j][i];
    for(j=0;j<ncl[i];j++) dumcctr[j][0]=qtl[j];
    kmeans(dumy,n,1,dumcctr,ncl[i],dum2,ning,maxiter,wss1,&kopt);
    for(j=0;j<ncl[i];j++) cctr[k++]=dumcctr[j][0];
    FREE_VECTOR(probs);
    FREE_VECTOR(wss1);
    FREE_VECTOR(ning);
    FREE_VECTOR(qtl);
    FREE_MATRIX(dumcctr);
  }
  FREE_VECTOR(dum2);
  FREE_MATRIX(dumy);
  FREE_VECTOR(dum1);

  MAKE_VECTOR(cum1,dmin);
  cum1[0]=0;
  for(i=1;i<dmin;i++) cum1[i]=cum1[i-1]+ncl[i-1];
  MAKE_MATRIX(ccent,totcl,dmin);

  MAKE_VECTOR(buf,dmin);
  for (i=0; i<totcl; i++) {
    break_down(i, ncl, buf, dmin);
    for(j=0;j<dmin;j++) ccent[i][j]=cctr[cum1[j]+buf[j]];
  }
  FREE_VECTOR(cum1);
  FREE_VECTOR(cctr);
  FREE_VECTOR(buf);
  FREE_VECTOR(ncl);

  cent=eliminulls(y,n,dmin,&totcl,ccent,0);

  FREE_MATRIX(ccent);

  /*  printf("dmin = %d SVD %d\n",dmin,totcl);*/
  if(totcl>=nclus) {
    ind=0;
    MAKE_VECTOR(counts,totcl);
    MAKE_VECTOR(wss1,totcl);
    kmeans(y,n,dmin,cent,totcl,grpids,counts,maxiter,wss1,&k);
    FREE_VECTOR(wss1);
    FREE_VECTOR(counts);
    
    MAKE_VECTOR(grclass,totcl);
    hclassify(totcl,dmin,cent,2,nclus,grclass);
    
    MAKE_MATRIX(mu,totcl,dmin);
    MAKE_MATRIX(ltsigma,totcl,dmin*(dmin+1)/2);
    i=initials(cent,totcl,dmin,nclus,ningrp,mu,ltsigma,grclass);
    FREE_VECTOR(grclass);
    FREE_MATRIX(cent);
    for(i=0;i<n;i++) grpids[i]=assign_closest(y[i],dmin,nclus,mu);
    FREE_MATRIX(mu);
    FREE_MATRIX(ltsigma);
  }
  else FREE_MATRIX(cent);
  return ind;
}  


int starts_in_original_domain(int n,int m,double **x,int nclus,int *ningrp,
			      int *grpids)
{    
  double **cent,*dum1,*probs,*qtl,*cctr,**dumy,**dumcctr,
    *wss1,**ccent,**mu,**ltsigma,**y;
  int i,j,k,*counts,*ncl,totcl=1,sumcl=0,*dum2,*buf,*ning,kopt,
    *grclass,*cum1,maxiter=100000,ind=1;   

  MAKE_VECTOR(ncl,m);
  i=(int)ceil(pow(nclus,1.0/m));
  for(j=0;j<m;j++) ncl[j]=i;
  for(j=0;j<m;j++) {
    totcl*=ncl[j];
    sumcl+=ncl[j];
  }

  MAKE_MATRIX(y,n,m);
  MAKE_VECTOR(cctr,m);
  MAKE_VECTOR(dum1,m*(m+1)/2);
  meandispersion(x,n,m,cctr,dum1);
  
  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) y[i][j]=(x[i][j]-cctr[j])/sqrt(dum1[j*(j+1)/2+j]);
  }
  FREE_VECTOR(dum1);
  FREE_VECTOR(cctr);

  MAKE_VECTOR(dum1,n);
  MAKE_VECTOR(cctr,sumcl);
  MAKE_MATRIX(dumy,n,1);
  MAKE_VECTOR(dum2,n);

  k=0;
  for(i=0;i<m;i++) {
    MAKE_VECTOR(ning,ncl[i]);
    MAKE_VECTOR(qtl,ncl[i]);
    MAKE_MATRIX(dumcctr,ncl[i],1);
    MAKE_VECTOR(wss1,ncl[i]);
    for(j=0;j<n;j++) dum1[j]=y[j][i];
    MAKE_VECTOR(probs,ncl[i]);
    for(j=0;j<ncl[i];j++) probs[j]=(2*j+1)/(2.*ncl[i]);
    j=quantile(n,dum1,probs,qtl,ncl[i]);
    for(j=0;j<n;j++) dumy[j][0]=y[j][i];
    for(j=0;j<ncl[i];j++) dumcctr[j][0]=qtl[j];
    kmeans(dumy,n,1,dumcctr,ncl[i],dum2,ning,maxiter,wss1,&kopt);
    for(j=0;j<ncl[i];j++) cctr[k++]=dumcctr[j][0];
    FREE_VECTOR(probs);
    FREE_VECTOR(wss1);
    FREE_VECTOR(ning);
    FREE_VECTOR(qtl);
    FREE_MATRIX(dumcctr);
  }
  FREE_VECTOR(dum2);
  FREE_MATRIX(dumy);
  FREE_VECTOR(dum1);

  MAKE_VECTOR(cum1,m);
  cum1[0]=0;
  for(i=1;i<m;i++) cum1[i]=cum1[i-1]+ncl[i-1];
  MAKE_MATRIX(ccent,totcl,m);

  MAKE_VECTOR(buf,m);
  for (i=0; i<totcl; i++) {
    break_down(i, ncl, buf, m);
    for(j=0;j<m;j++) ccent[i][j]=cctr[cum1[j]+buf[j]];
  }
  FREE_VECTOR(cum1);
  FREE_VECTOR(cctr);
  FREE_VECTOR(buf);
  FREE_VECTOR(ncl);

  cent=eliminulls(y,n,m,&totcl,ccent,0);
  FREE_MATRIX(ccent);

  if(totcl>=nclus) {
    ind=0;
    MAKE_VECTOR(counts,totcl);
    MAKE_VECTOR(wss1,totcl);
    kmeans(y,n,m,cent,totcl,grpids,counts,maxiter,wss1,&k);
    FREE_VECTOR(wss1);
    FREE_VECTOR(counts);
    
    MAKE_VECTOR(grclass,totcl);
    hclassify(totcl,m,cent,2,nclus,grclass);
    
    MAKE_MATRIX(mu,totcl,m);
    MAKE_MATRIX(ltsigma,totcl,m*(m+1)/2);
    i=initials(cent,totcl,m,nclus,ningrp,mu,ltsigma,grclass);
    FREE_VECTOR(grclass);
    FREE_MATRIX(cent);
    for(i=0;i<n;i++) grpids[i]=assign_closest(y[i],m,nclus,mu);
    FREE_MATRIX(mu);
    FREE_MATRIX(ltsigma);
  }
  else FREE_MATRIX(cent);
  FREE_MATRIX(y);
  return ind;
}

int starts(int n,int m,double **y,int nclus,int dmin,double *D,int *ningrp,
	   int *grpids)
{
  /* This is a wrapper to whether we call something in the SVD domain or in the
     original domain. 
     In-built function: no use outside this program. */

  if (dmin==0) return starts_in_original_domain(n,m,y,nclus,ningrp,grpids);
  else
    return starts_in_svd_domain(n,m,y,nclus,dmin,D,ningrp,grpids);
}


int starts_svd(int n,int m,double **Mu,double **x,int nclus,int *ningrp,
		    double *pi,int *grpids,double **LTSigma,double alpha,
		    int llhdnotW)
{
  /*This is the main function which goes through the many choices.
    Note: alpha is the proportion of variance explained by the first dmin 
    singular values. (we use 0.999)
    llhdnotW = 1 means that the log-likelihood will be minimized over the terms
    of the singular values.
    0 implies W will be computed at the start points
*/ 
  double **y,**LTSigma_protem,sum,sum1,*pi_protem,**U,*D,**V,**Mu_protem,tm0, 
    tm1=Inf;
  int i,j,dmin,dmaxmin,gind=1,ind;   

  MAKE_MATRIX(U,n,m);  /* we assume n > m */
  MAKE_VECTOR(D,m);
  MAKE_MATRIX(V,m,m);

  i=svdd(x,n,m,D,U,V);
  FREE_MATRIX(V);
  sum=0;
  for(i=0;i<m;i++) sum+=SQ(D[i]);

  sum1=0;
  for (i=0;((sum1<alpha) && (i<m));i++) sum1+=D[i]*D[i]/sum;

  dmaxmin=i;

  MAKE_MATRIX(Mu_protem,nclus,m);
  MAKE_MATRIX(LTSigma_protem,nclus,m*(m+1)/2);
  MAKE_VECTOR(pi_protem,nclus);
  
  for(dmin=0;dmin<=dmaxmin;dmin++) {
    ind=0;
    if ((dmin==0)  && (dmaxmin<m)) dmin++; 
                                 /*no checking in the original space*/
    if (dmin==0) {
      MAKE_MATRIX(y,n,m);
      for(i=0;i<n;i++) {
	for(j=0;j<m;j++) y[i][j]=x[i][j];
      }
    }
    else  {
      MAKE_MATRIX(y,n,dmin);
      for(i=0;i<n;i++) {
	for(j=0;j<dmin;j++) y[i][j]=U[i][j];
      }
    }
    if(!starts(n,m,y,nclus,dmin,D,ningrp,grpids)) {
      FREE_MATRIX(y);
      i=initials(x,n,m,nclus,ningrp,Mu_protem,LTSigma_protem,grpids);
      
      if (llhdnotW) {
	for(i=0;i<nclus;i++) {
	  if (ningrp[i]<=m) ind=1;
	  pi_protem[i]=ningrp[i]/(1.*n);
	}
	if(ind==1) {
	  double *temp;
	  int kk=0;
	  MAKE_VECTOR(temp,m*(m+1)/2);
	  for(i=0;i<m*(m+1)/2;i++) temp[i]=0.;
	  for(j=0;j<nclus;j++) {
	    if (ningrp[j]>m) {
	      for(i=0;i<m*(m+1)/2;i++) temp[i]+=LTSigma_protem[j][i];
	      kk++;
	      ind=0;
	    }
	  }
	  for(j=0;j<nclus;j++) {
	    if (ningrp[j]<=m) {
	      for(i=0;i<m*(m+1)/2;i++) LTSigma_protem[j][i]=temp[i]/kk;
	    }	  
	  }
	  FREE_VECTOR(temp);
	}
	tm0=-lnlikelihood(n,m,nclus,pi_protem,x,Mu_protem,LTSigma_protem);
	/*	printf(" negllhd = %f \n",tm0);*/
	if (tm0<tm1) {
	  cpy(Mu_protem,nclus,m,Mu);
	  cpy(LTSigma_protem,nclus,m*(m+1)/2,LTSigma);
	  for(i=0;i<nclus;i++) pi[i]=pi_protem[i];
	  tm1=tm0;
	}
	gind=0;
      }
      else {
	gind=0;
	tm0=detw(nclus,m,LTSigma_protem,ningrp);
	if (tm0<tm1) {
	  cpy(Mu_protem,nclus,m,Mu);
	  tm1=tm0;
	}
      }
    }
    else FREE_MATRIX(y);
  }
  FREE_MATRIX(Mu_protem);
  FREE_MATRIX(LTSigma_protem);
  FREE_VECTOR(pi_protem);
  FREE_MATRIX(U);
  FREE_VECTOR(D);
  return gind;
}

int starts_via_svd(int n,int m,double **Mu,double **x,int nclus,int *ningrp,
		    double *pi,int *grpids,double **LTSigma,double alpha,
		    int llhdnotW)
/* All this function does is center the data and call starts_svd(), and then 
   put back the centers to the means and the centers */
{
  int i,j,ind;
  double *mu,*ltsigma;
  MAKE_VECTOR(mu,m);
  MAKE_VECTOR(ltsigma,m*(m+1)/2);
  meandispersion(x,n,m,mu,ltsigma);
  FREE_VECTOR(ltsigma);
  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) x[i][j]-=mu[j];
  }
  ind=starts_svd(n,m,Mu,x,nclus,ningrp,pi,grpids,LTSigma,alpha,llhdnotW);
  if (!ind) {
    for(i=0;i<nclus;i++) {
      for(j=0;j<m;j++) Mu[i][j]+=mu[j];
    }
  }
  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) x[i][j]+=mu[j];
  }
  FREE_VECTOR(mu);
  return ind;
}
