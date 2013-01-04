/* mat_vec.c
 * 
 * Basic vector and matrix functions
 * 
 * Authors: Rouben Rostamian<rostamian@umbc.edu> and Ranjan Maitra<maitra@iastate.edu>
 * Fall 1996
 * Revised January 1999
 * Revised October 2000
 * Thoroughly revised March 2005
*/

#include <stdio.h>		/* for fprintf() */
#include <stdlib.h>		/* for malloc() */
#include <math.h>		/* for sqrt() */
#include <assert.h>
#include "array.h"

#include <R.h>
#include <Rinternals.h>

/* Allocates memory and creates an mxn Hilbert matrix of double */
double **dhilbert(int m, int n)
{
	int i, j;
	double **h;

	MAKE_MATRIX(h,m,n);
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			h[i][j] = 1.0 / ( 1.0 + i + j);

	return h;
}

/* Multiplies matrices a and b and puts the result in c which should be
 * pre-allocated.   exit() will be called if a and b are incompatible
*/
int multiply(double **a, int arows, int acols,
	      double **b, int brows, int bcols, double **c)
{
  int i, j, k;
  
  assert(acols==brows);
  
  for (i=0; i<arows; i++)
    for (j=0; j<bcols; j++) {
      c[i][j] = 0;
      for (k=0; k<acols; k++)
	c[i][j] += a[i][k] * b[k][j];
    }
  return 0;
}

/* Multiplies matrix a and vector x and puts the result in y which should be
 * pre-allocated.   exit() will be called if a and x are incompatible
*/
int matxvec(double **a, int arows, int acols,
		double *x, int xrows, double *y)
{
  int i, k;
  
  assert(acols==xrows);
  
  for (i=0; i<arows; i++){
    y[i] = 0;
    for (k=0; k<acols; k++){
      y[i] += a[i][k] * x[k];
    }
  }
  return 0;
}

/* Prints matrix with a spefified format */
void print_dmatrix(double **a, int rows, int cols, const char *format)
{
	int i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++)
			Rprintf(format, a[i][j]);
		//WCC putchar('\n');
		Rprintf("\n");
	}
}

/* Prints vector with a spefified format */
void print_dvector(double *a, int rows, const char *format)
{
  int i;
  for (i=0; i<rows; i++) {
    Rprintf(format, a[i]);
  }
  return;
}

void matrpose(double **a,int rows,int cols,double **aT)
{
  int i,j;
  for (i=0;i<rows;i++) {
    for (j=0;j<cols;j++) {
      aT[j][i]=a[i][j];
    }
  }
  return;
}

/* Calculates and returns the Euclidean norm of a vector.
 * Modeled after BLAS's norm2 routine.
 *
 * The `sum' variable accummulates (x_i/scale)^2, where `scale'
 * is the largest |x_i| seen so far.
*/
double dEnorm(double *x, int n)
{
	double scale = 0.0;
	double sum = 1.0;
	double b, r;
	int i;

	if (n < 1)
		return 0.0;
	else if (n == 1)
		return fabs(x[0]);

	for (i=0; i<n; i++) {

		if (x[i]==0.0)
			continue;

		b = fabs(x[i]);

		if (b < scale) {
			r = b/scale;
			sum += r*r;
		} else {
			r = scale/b;
			sum = sum*r*r + 1.0;
			scale = b;
		}
	}

	return scale * sqrt(sum);
}

/*copy matrix A to matrix B*/
int cpy(double **a,int nrows,int ncols,double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[i][j];
    }
  }
  return 0;
}

double quadratic(double **A,double *x,int p) 
{                  /*calculates x'Ax where A is a p x p-dimensional matrix */
  int i,j;
  double qform=0.0;
  for(i=0;i<p;i++) {
    for(j=0;j<p;j++) {
      qform+=x[i]*x[j]*A[i][j];
    }
  }
  return(qform);
}

double ltquadratic(double *ltA,double *x,int p) 
{                  /*calculates x'Ax where A is a p x p-dimensional symmetric 
		     matrix in packed lower-triangular form*/
  int i,j;
  double qform=0.0;
  for(i=0;i<p;i++) {
    qform+=x[i]*x[i]*ltA[((i+1)*(i+2))/2-1];
    for(j=0;j<i;j++) {
      qform+=2.*x[i]*x[j]*ltA[i*(i+1)/2+j];
    }
  }
  return(qform);
}


int ar(double **A,int n,double rho) /* sets up a AR-1 correlation matrix 
				       with rho */
{
  int i,j;
  double *x;
  MAKE_VECTOR(x,n);
  x[0]=1.0;
  for(i=1;i<n;i++) {
    x[i]=x[i-1]*rho;
  }
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      A[i][j]=x[abs(i-j)];
    }
  }
  FREE_VECTOR(x);
  return 0;
}

int arinv(double **A,int n,double rho) /* sets up the inverse of a AR-1
					  correlation matrix  with rho */
{
  int i,j;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      if (((i==0) && (j==0)) || ((i==(n-1)) && (j==(n-1)))) {
	A[i][i]=1./(1-rho*rho);
      }
      else {
	if (i==j) {
	  A[i][i]=(1+rho*rho)/(1-rho*rho);
	}
	else {
	  if (abs(i-j)==1) {
	    A[i][j]=-rho/(1-rho*rho);
	  }
	  else {
	    A[i][j]=0.;
	  }
	}
      }
    }
  }
  return 0;
}

