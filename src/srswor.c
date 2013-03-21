#include<stdlib.h>
#include "order.h"
#include "array.h"

/* Marked by Wei-Chen Chen on 2008/10/29
#define MATHLIB_STANDALONE 1
*/
                              /*It is essential to have this before the call 
                               to the Rmath's header file because this decides
                               the definitions to be set. */


#include <R.h>
#include <Rmath.h>

/* Note that use of this function involves a prior call to the Rmath library to
   get the seeds in place. It is assumed that this is done in the calling 
   function. */

/* Equal probability sampling; without-replacement case */
/* Adapted from the R function called SampleNoReplace */


int srswor(int n, int k, int *y)
{
	/* Provide k out of n indices sampled at random without replacement */
	
	if (k > n) {

		//WCC printf("Error: k = %d  greater than n = %d  in srswor()\n", k, n);
		REprintf("Error: k = %d  greater than n = %d  in srswor()\n", k, n);
		return 1;
	}
	else {
		
		int i, j;
		int *x;
		
		MAKE_VECTOR(x, n);
		for (i = 0; i < n; i++)	x[i] = i;
		
/* Add by Wei-Chen Chen on 2008/10/29
   This is required to get the initial seed state.  */
GetRNGstate();

		for (i = 0; i < k; i++) {
			j = (int) n * runif(0, 1);
			y[i] = x[j];
			x[j] = x[--n];
		}

/* Add by Wei-Chen Chen on 2008/10/29
   This is required to set the final seed state.  */
PutRNGstate();

		FREE_VECTOR(x);
	}
	return 0;
}


/* Not implemented function, marked by Wei-Chen Chen on 2008/09/28.
int WRSampleUnequalProb(int n, int k, double *prob, int *smpl)
{
*/
  /*
    Given n and k, provide a random sample with replacement of size n from 
    k classes, each occurring with unequal probabilities, or frequencies.
    Input parameters are as follows:

    n = sample size
    k = numer of class ids
    prob = probability of occurrence of class ids (also can handle frequencies)
    smpl = n-variate vector containing ids (0 through k-1) of the sample

    written by Ranjan Maitra, January 21, 2007.  
*/


/*
}
*/
