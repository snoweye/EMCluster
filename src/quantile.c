/**
 * Compute a given set of quantiles for a one-dimensional arrays
 *
 * Ranjan Maitra, Ames, IA 50014
 * 2005/09/15
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include"array.h"

#define RETURN_CMP(a,b) if (a < b) { return -1; } \
  else if (a == b) { return 0; } \
  else { return 1; }

int CompareDouble(const void* v1, const void* v2)
{
  /* written by David Faden. All rights reserved. */
  /* same as in sorted.c -- just want to write in contained form */

  double d1 = *(const double*)v1;
  double d2 = *(const double*)v2;
  RETURN_CMP(d1, d2);
}

double* copyArray(const double* values,int len)
{
  /* modified from program written by David Faden. All rights reserved. */
  double* result = malloc(sizeof(double)*len);
  memcpy(result, values, len * sizeof(double));
  return result;
}

int quantile(int n,double *x,double *p,double *q, int numqs)
/* Similar function to the Splus/R quantile function. Returns in q[i]
   the p[i]th quantile of the dataset in the n-dimensional vector x. Here, 
   numqs represent the number of    quantiles desired to be returned */
{
  double* sorted,resid;
  int i,modul;
  sorted = copyArray(x,n);

  qsort(sorted,n,sizeof(double),CompareDouble);
  
  for (i=0;i<numqs;i++) {
    modul=(int)floor(p[i]*(n-1));
    resid=p[i]*(n-1)-modul;
    if (resid==0)    q[i]=sorted[modul];
    else  q[i]=resid*sorted[modul+1]+(1-resid)*sorted[modul];
  }
  free(sorted);
  return 0;
}

double trimmed_mean(int n,double *x,double left,double right)
/* This function calculates the trimmed mean of x after dropping the lowest 
   "left" proportion of observations and the upper "right" proportion of 
   observations. The sum of "left" and "right, each assumed to be between 0 
   and 1, is further assumed to be between 0 and 1 and the code does not check
   for consistency.*/
{
  int i;
  double* sorted,nx=0;
  double sum = 0.0;

  sorted = copyArray(x,n);

  qsort(sorted,n,sizeof(double),CompareDouble);
  
  for (i=(int)left*n;i<(int)right*n;i++) {
    sum+=sorted[i];
    nx++;
  }
  free(sorted);
  return(sum/nx);
}

