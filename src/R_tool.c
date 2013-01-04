/* This file contains some tool functions required by other functions
   R_*() in "src/R_*.c".

   Writen: Wei-Chen Chen on 2008/09/28.
*/

#include<R.h>
#include<Rinternals.h>

/* Allocate a pointer array with double precision. */
double** allocate_double_array(int n){
  double **pointerarray = (double **) malloc(n * sizeof(double *));
  if(pointerarray == NULL){
    error("Memory allocation fails!\n");
  }
  return(pointerarray);
} /* End of allocate_double_array(). */
