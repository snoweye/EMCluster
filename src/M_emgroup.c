/* This file contains functions modified from "emgroup.c" using in
   Dr. Maitra's ten-clusters programs.

   Modified: Wei-Chen Chen on 2008/12/03.
*/

#include "mb_tool.h"

/* Main function for emgroup(). See "src/R_M_emgroup.c" for details */
int M_emgroup(double **x,int n,int p,int nclass,double *pi,double **Mu,
    double **LTSigma,double *llhdval,int *nc,int *class,
    double alpha, int em_iter, double em_eps,
    int *conv_iter, double *conv_eps){
  int j,flag=0;
  double like=0.0;
  
  if (nclass==1) {
    nc[0]=n;
    pi[0]=1.0;
    for (j=0;j<n;j++) class[j]=0;
/* These formulae are not correct with two errors,
   1. n-1 should be replace by n in meandispersion() for LTSigma, and
   2. determinant() will change the values of LTSigma[0], see "initials.c".
   Modified: Wei-Chen Chen on 2008/12/05.
    meandispersion(x,n,p,Mu[0],LTSigma[0]);
*/
    meandispersion_MLE(x,n,p,Mu[0],LTSigma[0]);
    like=-0.5*n*p-0.5*n*log(determinant(LTSigma[0],p))-0.5*n*p*log(2*M_PI);
  }
  else {
    if(!starts_via_svd(n,p,Mu,x,nclass,nc,pi,class,LTSigma,alpha,1))      {
      for(j=0;j<nclass;j++) pi[j]=nc[j]/(double)n;
      emcluster(n,p,nclass,pi,x,Mu,LTSigma,em_iter,em_eps,&like,
                conv_iter,conv_eps);
      assign(n,p,nclass,x,pi,Mu,LTSigma,class,nc);
    }
    else flag=1;
  }
  (*llhdval)=like;
  return flag;
} /* END of M_emgroup(). */

