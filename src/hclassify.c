/* This code is a wrapper to the hc.c and hclass.c code. What it does is 
   return the class indicators for a hierarchically clustered tree, using some 
   user-specified criterion (out of a list of seven) and the number of classes 
   nclass */

#include<stdio.h>
#include "array.h"

void hc(int n, int m, int iopt, double **data, int *ia, int *ib, double *crit);

void hclass(int n, int *ia, int *ib, int lev, int *iclass);

void hclassify(int n,int m, double **x,int hcrit,int nclass,int *class)
{
  double *crit;
  int *ia,*ib;
  MAKE_VECTOR(ia,n);
  MAKE_VECTOR(ib,n);
  MAKE_VECTOR(crit,n);

  hc(n,m,hcrit,x,ia,ib,crit);
  FREE_VECTOR(crit);
  
  hclass(n,ia,ib,nclass,class);

  FREE_VECTOR(ia);
  FREE_VECTOR(ib);
  return;
}
