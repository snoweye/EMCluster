#include<stdlib.h>
#include "array.h"

int sort(int n,double *x);

void unique(int n,double *x,int *m,double *y) /* Returns the unique values in 
						 x in ascending order: the 
						 number of unique values is 
						 output in m and the unique 
						 values are output in sorted
						 order in y */
{
  int i,j;
  double *dum;

  MAKE_VECTOR(dum,n);
  for(i=0;i<n;i++) dum[i]=x[i];
  sort(n,dum);
  i=0;
  j=0;
  while (j<n) {
  y[i]=dum[j];
  while ((j<n) && (dum[j]==y[i])) {
    j++;
  }
  i++;
  }
  (*m)=i;
  FREE_VECTOR(dum);
  return;
}

