#include <stdio.h> 
#include "array.h"

#include <R.h>
#include <Rinternals.h>

void dsyevr_(char *jobz, char *range, char *uplo,
	     int *n, double *a, int *lda,
	     double *vl, double *vu, int *il, int *iu,
	     double *abstol, int *m, double *w,
	     double *z, int *ldz, int *isuppz,
	     double *work, int *lwork, int *iwork, int *liwork,
	     int *info);

void dspevd_(char *jobz,char *UPLO,int *n,double *ap,double *w,double *z,
	     int *ldz,double *work,int *lwork,int *iwork,int *liwork,int *info);


/* C front end to LAPACK's dspevd_() routine.
 *
 * Calculates the eigenvalues and the eigenvectors of the nxn symmetric matrix 
 * A and returns them in the vector E, and the vector EV, respectively.
 * The i'the eigenvector corresponds to the i'th eigenvalue.
 *
 * Note that the eigenvalues are returned in the descending order.
 * 
 * */

int eigend(double *A, double *EV, double *E, int n)
{
  int lwork=-1,liwork=-1;
  double dx,*z,*a,*w,*work = &dx;/* work points to a temporary cell */
  int info,ix,i,j,*iwork = &ix;   /* iwork points to a temporary cell */
  char jobz='V', uplo='U';

  MAKE_VECTOR(a,n*(n+1)/2);
  MAKE_VECTOR(w,n);
  MAKE_VECTOR(z,n*n);

  for(i=0;i<n*(n+1)/2;i++) a[i]=A[i];

  /* Call dspevd_() with lwork=-1 and liwork=-1 to query the optimal sizes of 
   * the work and iwork arrays.
   * */
  dspevd_(&jobz,&uplo,&n,a,w,z,&n,work,&lwork,iwork,&liwork,&info);
  
  if (info==0) {
    lwork = (int)*work;
    liwork = *iwork;
    
    /* allocate optimal sizes for work and iwork */
    MAKE_VECTOR(work,lwork);
    MAKE_VECTOR(iwork,liwork);
    
    if (work!=NULL && iwork!=NULL) {
      dspevd_(&jobz,&uplo,&n,a,w,z,&n,work,&lwork,iwork,&liwork,&info);
      if (info==0) {
	for(i=0;i<n;i++) {
	  E[i]=w[n-1-i];
	  for (j=0;j<n;j++) {
	    EV[j*n+i]=z[(n-1-j)*n+i];
	  }
	}
      }
      else {
	REprintf("error in dspvd at calculation stage: Error code %d\n",info);
      }
    }
    FREE_VECTOR(work);
    FREE_VECTOR(iwork);
  }
  FREE_VECTOR(a);
  FREE_VECTOR(w);
  FREE_VECTOR(z);
  return info;
}

/* C front end to LINPACK's dsyevr_() routine.
 *
 * Calculates the eigenvalues and eigenvectors of the nxn symmetric matrix A.
 * The eigenvalues are returned in the vector w.
 * The (orthonormal) eigenvectors are returned in the matrix z.
 * The ith column of z holds the eigenvector associated with w[i].
 * Written by Rouben Rostamian and Ranjan Maitra
 * */
int LP_sym_eigvecs(double *a, int n, double *w, double *z)
{
        int lwork=-1, liwork=-1;
        double dx, *work = &dx;         /* work points to a temporary cell */
        int ix, *iwork = &ix;           /* iwork points to a temporary cell */
        double abstol = -1.0;           /* force default */
        char jobz='V', range='A', uplo='L';
        int il, iu;
        double vl, vu;
        int info, m;

        /* Call dsyevr_() with lwork=-1 and liwork=-1 to query the
         * optimal sizes of the work and iwork arrays.
         * */
        dsyevr_(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &iu,
                        &abstol, &m, w, NULL, &n, NULL, 
                        work, &lwork, iwork, &liwork, &info);

        if (info!=0) {
//WCC                fprintf(stderr, "trouble in file __FILE__, line __LINE__: "
                REprintf("trouble in file __FILE__, line __LINE__: "
                                "info = %d\n", info);
        }
	else {
	  int *isuppz;
	  lwork = (int)*work;
	  liwork = *iwork;
	  
	  /* allocate optimal sizes for work and iwork */
	  work = malloc(lwork*sizeof *work);
	  iwork = malloc(liwork*sizeof *iwork);
	  isuppz = malloc(2*n*sizeof *isuppz);
	  
	  if (work==NULL || iwork==NULL || isuppz==NULL) {
//WCC	    fprintf(stderr, "trouble in file __FILE__, line __LINE__: "
	    REprintf("trouble in file __FILE__, line __LINE__: "
		    "out of memory!\n");
	  }
	  else {
	    /* now call dsyevr_() in earnest */
	    dsyevr_(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &iu,
		    &abstol, &m, w, z, &n, isuppz, 
		    work, &lwork, iwork, &liwork, &info);
	    
	    if (info!=0) {
//WCC	      fprintf(stderr, "trouble in file __FILE__, line __LINE__: "
	      REprintf("trouble in file __FILE__, line __LINE__: "
		      "info = %d\n", info);
	    }
	  }
	  free(isuppz);
	}
        free(work);
        free(iwork);

        return info;
}

/* Front end to LP_sym_eigvecs
 *
 * Calculates the eigenvalues and eigenvectors of the nxn symmetric matrix A.
 * The eigenvalues are returned in the vector w.
 * The (orthonormal) eigenvectors are returned in the matrix z.
 * The ith column of z holds the eigenvector associated with w[i].
 * Written by Ranjan Maitra
 * Note that one major task done here is to reverse the order of the 
 * eigenvalues (to something that makes more sense) and to put them in
 * decreasing order and the corresponding eigenvectors
 * */
int symeigens(double *a, int n, double *w, double *z)
{
  double *eval,*evec;
  int  i,j,info;
  MAKE_VECTOR(eval,n);
  MAKE_VECTOR(evec,n*n);
  info=LP_sym_eigvecs(a,n,eval,evec);
  if (info==0) {
    for(i=0;i<n;i++) {
      w[i]=eval[n-1-i];
      for (j=0;j<n;j++) {
	z[j*n+i]=evec[(n-1-j)*n+i];
      }
    }
  }
  FREE_VECTOR(eval);
  FREE_VECTOR(evec);
  return info;
}

/*
  Calculates eigevectors and eigenvalues for a packed symmetric matrix. 
  The eigenvalues returned are in descending order.
  The returned eigenvectors correspond to the eigenvalues.
 */

int eigens(double *A, double *EVec, double *EVal, int n)
{
  double *a;
  int i,j;
  MAKE_VECTOR(a,n*n);
  for(i=0;i<n;i++) {
    for(j=0;j<i;j++) {
      a[i*n+j]=A[i*(i+1)/2+j];
      a[j*n+i]=A[i*(i+1)/2+j];
    }
    a[i*n+i]=A[i*(i+1)/2+i];
  }
  i=symeigens(a,n,EVal,EVec);
  FREE_VECTOR(a);
  return i;
}
