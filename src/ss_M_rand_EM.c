/* This file is modified from "M_rand_EM.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"


/* This function is coded in "src/ss_M_rand_EM.c".
   Input:
     **x: double[n, p], data matrix of n*p.
     n: int[1], number of observations.
     p: int[1], number of dimersions.
     nclass: int[1], number of classes.
     *pi: double[nclus], proportions of clusters.
     **Mu: double[nclass, p], means of MVNs.
     **LTSigma: double[nclus, p * (p + 1) / 2], lower triangular sigma matrices.
     *llhdval: array[1], log likelihood value.
     *nc: int[nclass], number of observations in each class.
     shortiter: int[1], number of short iterations, 500 by default.
     fixediter: int[1], number of rand iterations, 1 by default.
     em_iter: int[1], number of emclust() iterations, 1000 by default.
     em_eps: double[1], epsilon for emclust() iterations, 1e-4 by default.
     conv_iter: int[1], number of convertent iterations.
     conv_eps: double[1], epsilon at convertent iterations.
     lab: int[n], -1 for points with unknown clusters; 0,..,(labK-1) for known.
     labK: int[1], the number of known clusters.
   Output:
     flag: int[1], a returned flag from rand_EM().
     *pi, **Mu, **LTSigma, *llhdval, *nc, *class, *conv_iter, *conv_eps.
*/
int ss_M_rand_EM(double **x, int n, int p, int nclass, double *pi,
    double **Mu, double **LTSigma, double *llhdval, int *nc, int *class,
    int shortiter, int fixediter,
    int em_iter, double em_eps, int *conv_iter, double *conv_eps,
    int *lab, int labK){
  int j;

  if(nclass == 1){
    nc[0] = n;
    pi[0] = 1.0;
    for(j = 0; j < n; j++) class[j] = 0;
    meandispersion_MLE(x, n, p, Mu[0], LTSigma[0]);
    *llhdval = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
               0.5 * n * p * log(2 * PI);
  } else{
    ss_mod_shortems(n, p, nclass, pi, x, Mu, LTSigma,
                    shortiter, fixediter, conv_iter, conv_eps,
                    lab, labK);
    ss_emcluster(n, p, nclass, pi, x, Mu, LTSigma, em_iter, em_eps, llhdval,
                 conv_iter, conv_eps, lab);
    ss_assign(n, p, nclass, x, pi, Mu, LTSigma, class, nc, lab);
  } 

  return 0;
} /* End of ss_M_rand_EM(). */

