/* Matrix inverse. Use posymatinv for inversion of a positive-definite 
   matrix. 
  
   Author: Ranjan Maitra <maitra@iastate.edu>
   Date:   03/07/2005
   Uses:   LAPACK 
*/  

#include <stdio.h>              /* for fprintf() */
#include <stdlib.h>
#include <string.h>
#include "array.h"
#include "mat_vec.h"

void dgetrf_(int *Mp, int *Np, double *A, int *LDA, int *PIVOT, int *INFOp);
void dgetri_(int *Np, double *A, int *LDA, int *PIVOT, double *WORK, 
	     int *LWORK, int *INFOp);

void  dpotrf_(char *UPLOp,int *Np, double *A, int *LDAp, int *INFOp);
void  dpotri_(char *UPLOp,int *Np, double *A, int *LDAp, int *INFOp);

void dpptrf_(char *UPLOp,int *Np,double *A,int *INFOp);
void dpptri_(char *UPLOp,int *Np,double *A,int *INFOp);

int matinv(int sizeA,double **A,double (*determinant))
{
  int i, j , *pivot,N=sizeA*sizeA,size=sizeA;
  double *AT,*work;	/* AT=transpose vectorized matrix (to accomodate
			   Fortran) 
			   work=workspace vector */
  int INFO,ipiv=1;

  MAKE_VECTOR(AT,size*size);
  MAKE_VECTOR(work,size*size);
  MAKE_VECTOR(pivot,size);

  for (i=0; i<size; i++)		/* to call a Fortran routine from C */
    {				/* have to transform the matrix */
      for(j=0; j<size; j++) AT[j+size*i]=A[j][i];		
    }						
  
  dgetrf_(&size,&size,AT,&size,pivot,&INFO);
  /* LAPACK routine DGETRF computes an LU factorization of a general 
     m x n matrix A using partial pivoting with row interchanges. The
     factorization has the form A = P * L * U where P is a permutation
     matrix, L is lower triangular with unit diagonal elements (lower 
     trapezoidal if m > n), and U is upper triangular (upper trapezoidal
     if m < n). Note that because of the permutation, the determinant 
     needs to be multiplied by -1 for every interchange that has occurred.
     Parameters in the order as they appear in the function call: 
     number of rows of the matrix A, number of columns of the 
     matrix A, the matrix A, the leading dimension of A, the 
     array that records pivoting, and the flag for the
     result. On exit, A contains the factors of L and U (with the
     diagonals of L not stored).*/	  
  if (INFO==0) {
    for(i=0;i<size;i++) {
      if (i!=(pivot[i]-1)) ipiv*=-1; /* PIVOT assumes indices are from 1 
					through N*/
    }
    (*determinant)=(double)ipiv;
    for (i=0;i<size;i++) {
      (*determinant)*=AT[i+i*size];
    }
    dgetri_(&size,AT,&size,pivot,work,&N,&INFO); 
    /* LAPACK routine DGETRI computes the inverse of a matrix A 
       using the output of DGETRF. This method inverts U and then
       computes A^(-1) by solving A^(-1)L = U^(-1) for A^(-1).
       parameters in the order as they appear in the function call:
       order of the matrix A, the matrix A, the leading dimension of
       A, the array that records pivoting, workspace, the 
       dimension of the workspace array, and the flag for the 
       result. On exit, A contains the inverted matrix. */
    if (INFO!=0) {
/* Marked by Wei-Chen Chen on 2009/06/07.
*     printf("Problem in matinv: dgetri error %d\n",INFO);
*/
    }
  }
  else {
/* Marked by Wei-Chen Chen on 2009/06/07.
*   printf("Problem in matinv: dgetrf error %d\n",INFO);
*/
  }
  for (i=0; i<size; i++)		/* to call a Fortran routine from C */
    {				/* have to transform the matrix */
      for(j=0; j<size; j++) {
	A[j][i]=AT[j+size*i];		
      }
    }
  FREE_VECTOR(AT);
  FREE_VECTOR(pivot);
  FREE_VECTOR(work);
  return 0;
}

int posymatinv(int size,double **A,double (*determinant))
{
  int i, j,INFO,N,LDA ;
  char uplo='L';
  double *AT;  /* AT=transpose vectorized matrix (to accomodate Fortran) */
  
  MAKE_VECTOR(AT,size*size);
  for (i=0; i<size; i++)		/* to call a Fortran routine from C */
    {				/* have to transform the matrix */
      for(j=0; j<size; j++) AT[j+size*i]=A[j][i];		
    }						
  
  N=size;
  LDA=size;

  dpotrf_ (&uplo, &N, AT, &LDA, &INFO);
  /* LAPACK routine DPOTRF computes an Cholesky decomposition of
     a symmetric positive definite matrix A. 
     Parameters in the order as they appear in the function call:
     uplo="U" indicates that the strictly lower triangular part of
     A will be ignored, N is the order of the matrix A, the 
     matrix A, the leading dimension of A, and the flag for the
     result. On exit, the upper triangle of A contains U.*/  
  if (INFO==0) {
    int i;
    (*determinant)=1.0;
    for (i=0;i<N;i++) {
      (*determinant)*=AT[i+i*N]*AT[i+i*N];
    }
    dpotri_ (&uplo, &N, AT, &LDA, &INFO);
    /* LAPACK routine DPOTRI computes the inverse of a matrix A 
       using the output of DPOTRF. This method inverts U using the 
       Cholesky factorization of A. 
       Parameters in the order as they appear in the function call:
       uplo="U" indicates that the strictly lower triangular part of
       A will be ignored, c1 is the order of the matrix A, the 
       matrix A, the leading dimension of A, and the flag for the
       result. On exit, the upper triangle of A contains U.*/  
    if (INFO!=0) {
/* Marked by Wei-Chen Chen on 2009/06/07.
*     printf("Problem in posymatinv: dpotri error %d\n",INFO);
*/
    }
  }
  else {
/* Marked by Wei-Chen Chen on 2009/06/07.
*   printf("Problem in posymatinv: dpotrf error %d\n",INFO);
*/
  }

  for (i=0; i<size; i++) {    /*transform the matrix back*/
      for(j=i; j<size; j++) {
	A[j][i]=AT[j+size*i];		
	A[i][j]=AT[j+size*i];		
      }
    }						
  FREE_VECTOR(AT);
  return 0;
}

int pposymatinv(int N,double *A, char UPLO, double *determinant)
{
  int INFO;

  dpptrf_(&UPLO,&N,A,&INFO);
  /* LAPACK routine DPPTRF computes the Cholesky factorization of 
     a packed real symmetric positive definite matrix A stored in 
     packed format. The factorization has the form 
     A = U**T * U,  if UPLO = 'U', or  A = L  * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular.
     Parameters in the order in which they appear in the function 
     call:
     uplo="U" indicates the upper triangle, while uplo="L" indicates 
     the lower triangle of the matrix packed into the vector A, N is 
     the order of the matrix, A is the vector containing the packed 
     matrix and INFO is the flag for the result. On exit, if INFO = 0,
     the triangular factor U or L from the Cholesky factorization 
     A = U**T*U or A = L*L**T, in the same storage format as A. */
  if (INFO==0) {
    int i;
    (*determinant)=1.0;
    if (UPLO=='U') {
      for (i=0;i<N;i++) {
        (*determinant)*=A[i+i*(i+1)/2]*A[i+i*(i+1)/2];
      }
    }
    else {
      for (i=0;i<N;i++) {
        (*determinant)*=A[i+(i*(2*N-i-1))/2]*A[i+(i*(2*N-i-1))/2];
      }
    }
    dpptri_(&UPLO,&N,A,&INFO);
    /* LAPACK routine DPPTRI computes the inverse of a packed real 
       symmetric positive definite matrix A using the Cholesky 
       factorization A = U**T*U or A = L*L**T computed by DPPTRF. 
       Parameters in the order in which they appear in the function 
       call:
       uplo="U" indicates the upper triangle, while uplo="L" indicates
       the lower triangle of the matrix packed into the vector A, N is 
       the order of the matrix, A is the vector containing the packed 
       matrix and INFO is the flag for the result. On exit, if INFO = 0, 
       the upper or lower triangle of the (symmetric) inverse of A, 
       overwriting the input factor U or L.*/
    if (INFO!=0) {
/* Marked by Wei-Chen Chen on 2008/12/26.
*     printf("\rProblem in pposymatinv: dpptri error %d",INFO);
*     Rprintf("dpptri error %d\n",INFO);
*/
    }
  }
  else {
/* Marked by Wei-Chen Chen on 2008/12/26.
*   printf("\rProblem in pposymatinv: dpptrf error %d",INFO);
*   Rprintf("dpptrf error %d\n", INFO);
*/
  }
  return  INFO;
}

