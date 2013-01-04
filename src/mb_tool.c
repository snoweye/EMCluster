/* This file contains several functions to perform model-based initializers.
   Created: Wei-Chen Chen on 2009/03/04.
*/

#include "mb_tool.h"

double mb_median(int n, double *x){
  double p = 0.5, q = 0.0;
  quantile(n, x, &p, &q, 1);
  return(q);
} /* End of mb_median(). */
double mb_quantile(int n, double *x, double p){
  double q = 0.0;
  quantile(n, x, &p, &q, 1);
  return(q);
} /* End of mb_quantile(). */


double fcn(double lambda, void *pt_param){
  PARAM *param = pt_param;
  return ppois(param->n_nclass, lambda, 1, 0) - 1.0 + 1e-6;
} /* End of fcn(). */
double find_lambda(void *pt_param){
  double lambda = 0.0;
  PARAM *param = pt_param;

  lambda = findzero(param->lower_bound, param->upper_bound, fcn, pt_param,
                    &param->tol, &param->maxit);

  if(param->maxit == -1){
    //WCC printf("find_lambda() fails.\n"); 
    //WCC exit(1);
    error("find_lambda() fails.\n");
  }

  return lambda;
} /* End of find_lambda(). */

