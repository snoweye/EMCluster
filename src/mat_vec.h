#ifndef MATVEC_H
#define MATVEC_H

double **identmatrix(int n);
double dEnorm(double v[], int n);
int multiply(double **a, int arows, int acols,
		double **b, int brows, int bcols, double **c);
void print_dmatrix(double **a, int rows, int cols, const char *format);

void print_dvector(double *a, int rows, const char *format);
int matxvec(double **a, int arows, int acols,
	    double *x, int xrows, double *y);
void matrpose(double **a,int rows,int cols,double **aT);
int arinv(double **A,int n,double rho); /* sets up the inverse of a AR-1
					  correlation matrix  with rho */
int ar(double **A,int n,double rho); /* sets up a AR-1 correlation matrix 
				       with rho */
double quadratic(double **A,double *x,int p); 
     /*calculates x'Ax where A is a p x p-dimensional matrix */
double ltquadratic(double *A,double *x,int p); 
     /*calculates x'Ax where A is a p x p-dimensional symmetric matrix in 
      packed lower-triangular form*/
int cpy(double **a,int nrows,int ncols,double **b);

#endif /* MATVEC_H */
