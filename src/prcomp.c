#include "array.h"
#include "math.h"
#include "mat_vec.h"

int svdd(double **a, int m, int n, double *d, double **u, double **v);
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);

void prcomp(int n, int m,double **x,double **UtX,double *D)
{
  int i,j;
  double *mu,*nu,*ltsigma,**xx,**V,**Vt;

  MAKE_VECTOR(mu,m);
  MAKE_VECTOR(ltsigma,m*(m+1)/2);  
  meandispersion(x,n,m,mu,ltsigma);
  FREE_VECTOR(ltsigma);
  MAKE_MATRIX(xx,n,m);
  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) xx[i][j]=(x[i][j]-mu[j])/sqrt(n-1.);
  }

  MAKE_MATRIX(V,m,m);
  i=svdd(xx,n,m,D,UtX,V);
  MAKE_MATRIX(Vt,m,m);
  matrpose(V,m,m,Vt);
  MAKE_VECTOR(nu,m);
  i=matxvec(Vt,m,m,mu,m,nu);
  FREE_MATRIX(Vt);
  FREE_VECTOR(mu);
  FREE_MATRIX(V);
  FREE_MATRIX(xx);

  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) {
      UtX[i][j]*=D[j];
      UtX[i][j]+=nu[j];
    }
  }
  FREE_VECTOR(nu);
  return;
}

