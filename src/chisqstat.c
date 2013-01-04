#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mat_vec.h"
#include "array.h"
#define Inf 1e+140

int posymatinv(int size,double **A,double (*determinant));
int pposymatinv(int N,double *A, char UPLO, double (*determinant));

double chisqstatinv(int p, double *X,double *mu,double **A)
{
  /* This function computes the test statistic (X-mu)'A(X-mu) also 
     given by the chisquare-statistic for the normal distribution. */
  int i;
  double *y,temp;
  MAKE_VECTOR(y,p);
  for (i=0;i<p;i++) {
    y[i]=X[i]-mu[i];
  }
  temp=quadratic(A,y,p);
  FREE_VECTOR(y);
  return(temp);
}

double chisqstat(int p, double *X,double *mu,double **Sigma,double *detSig)
{
  /* This function computes the test statistic (X-mu)'Sigma^{-1}(X-mu) also 
     given by the chisquare-statistic for the normal distribution. It also 
     returns the determinant of Sigma for use in other programs. All it does 
     is calculates the inverse of Sigma and the determinant and then transfer 
     control to chisqstatinv which is similar to the above but assumes that
     the inverse matrix and the determinant are supplied so that inverses 
     which can be easily calculated can be accomodated without going through
     the inversion step. */
  double **a,temp;
  int i;
  MAKE_MATRIX(a,p,p);
  cpy(Sigma,p,p,a);
  i=posymatinv(p,a,&(*detSig));
  temp=chisqstatinv(p,X,mu,a);
  FREE_MATRIX(a);
  return(temp);
}

double chisqstatltinv(int p, double *X,double *mu,double *ltA)
{
  /* This function computes the test statistic (X-mu)'A(X-mu) also 
     given by the chisquare-statistic for the normal distribution. Here ltA 
     contains the inverse of A in lower triangular form. */
  int i;
  double *z,temp;
  MAKE_VECTOR(z,p);
  for (i=0;i<p;i++) {
    z[i]=X[i]-mu[i];
  }
  temp=ltquadratic(ltA,z,p);
  FREE_VECTOR(z);
  return(temp);
}

double chisqstatlt(int p, double *X,double *mu,double *ltSigma,double *detSig)
{
  /* This function computes the test statistic (X-mu)'Sigma^{-1}(X-mu) also 
     given by the chisquare-statistic for the normal distribution. It also 
     returns the determinant of Sigma for use in other programs. It is assumed 
     that Sigma is symmetric and in packed lower triangular form. All it does 
     is calculates the inverse of Sigma and the determinant and then transfer 
     control to chisqstatltinv which is similar to the above but assumes that
     the inverse matrix and the determinant are supplied so that inverses 
     which can be easily calculated can be accomodated without going through
     the inversion step. */
  double *a,temp;
  int i;
  MAKE_VECTOR(a,p*(p+1)/2);
  for (i=0;i<p*(p+1)/2;i++) a[i]=ltSigma[i];
  i=pposymatinv(p,a,'U',&(*detSig));
  temp=chisqstatltinv(p,X,mu,a);
  FREE_VECTOR(a);
  return(temp);
}
