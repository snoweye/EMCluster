#include<stdlib.h>
#include "order.h"
#include "array.h"

#define RETURN_CMP(a,b) if (a < b) { return -1; } \
  else if (a == b) { return 0; } \
  else { return 1; }

int compareDouble(const void* v1, const void* v2)
{
  /* written by David Faden. All rights reserved. */

  double d1 = *(const double*)v1;
  double d2 = *(const double*)v2;
  RETURN_CMP(d1, d2);
}

int sort(int n,double *x)
/* This function sorts the n-dimensional array x in increasing order. It uses
   the standard library function qsort().
   The input array is replaced by the sorted array on output. 
   
   Written by Ranjan Maitra, Ames, IA 50014.
   2005/09/16. All rights reserved. */
{
  qsort(x,n,sizeof(double),compareDouble);
  return 0;
}

int mdimsort(int n,int p,double **x,int sortdim)
/* This function sorts the columns of the elements of an n x p matrix, in
   increasing order in the sortdim'th dimension. The input array is replaced 
   by the sorted array on output. 

   Written by Ranjan Maitra, Ames, IA 50014.
   2005/09/16. All rights reserved. */
   
{
  double *y,**xx;
  int i,j;
  size_t *ord;

  MAKE_VECTOR(y,n);
  for(i=0;i<p;i++) y[i]=x[i][sortdim];

  MAKE_VECTOR(ord,p);
  ord=orderDouble(y,n);

  MAKE_MATRIX(xx,n,p);
  for(i=0;i<n;i++) {
    for (j=0;j<p;j++) xx[i][j]=x[ord[i]][j];
  }
  for(i=0;i<n;i++) {
    for (j=0;j<p;j++) x[i][j]=xx[i][j];
  }
  FREE_MATRIX(xx);
  FREE_VECTOR(y);
  FREE_VECTOR(ord);
  return 0;
}
