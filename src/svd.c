#include <stdio.h> 
#include <stdlib.h> 
#include "array.h"

#define MIN(i,j) ((i)<(j) ? (i) : (j))
#define MAX(i,j) ((i)>(j) ? (i) : (j))

#include <R.h>
#include <Rinternals.h>

void dgesdd_(char *jobz,int *m,int *n,double *a,int *lda,double *s,double *u,
	     int *ldu,double *vt,int *ldvt,double *work,int *lwork,int *iwork,
	     int *info);
int svdd(double **a, int m, int n, double *d, double **u, double **v)
{
	double *A, *U, *VT;
	int lwork = -1;
	int liwork = 8*MIN(m,n);
	char jobz = 'S';
	double dw,*work=&dw; /*points to a temporary cell*/
	int *iwork;
	int i, j, k, info,minmn=MIN(m,n);

	MAKE_VECTOR(A, m*n);
	for (j=0, k=0; j<n; j++) {
		for (i=0; i<m; i++) A[k++] = a[i][j];
	}

	MAKE_VECTOR(U, m*minmn);
	MAKE_VECTOR(VT,minmn*n);
	MAKE_VECTOR(iwork, liwork);

	lwork=-1;
	
	dgesdd_(&jobz, &m, &n, A, &m, d, U, &m, VT, &n,
		work, &lwork, iwork, &info); /*call to get optimal lwork*/
	
	if (info!=0) {
	  //WCC printf("error: allocating LWORK in svdd\n");
	  //WCC exit(1);
	  error("error: allocating LWORK in svdd\n");
	}

	lwork=(int)*work;
	MAKE_VECTOR(work, lwork);

	dgesdd_(&jobz, &m, &n, A, &m, d, U, &m, VT, &minmn,
			work, &lwork, iwork, &info);
	FREE_VECTOR(A);
	FREE_VECTOR(work);
	FREE_VECTOR(iwork);

	for (j=0, k=0; j<minmn; j++) {
	  for (i=0; i<m; i++) u[i][j]=U[k++];
	}

	/* VT, as calculated by dgesdd_(), is the transpose of the right
	 * multiplier.  Here we undo the transpose so that the matrix
	 * v[][] returned by this function is not transposed anymore.
	 */

	for (i=0,k=0; i<n; i++) {
	  for (j=0; j<minmn; j++)   v[i][j]=VT[k++];
	}

	FREE_VECTOR(U);
	FREE_VECTOR(VT);

	return info;
}


