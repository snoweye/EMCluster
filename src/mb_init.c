/* This file contains several functions to perform model-based initializers.
   Created: Wei-Chen Chen on 2009/02/07.
*/

#include "mb_tool.h"

void cut_sub(double **X, int n, int p, int G, int min_n, double lambda,
    double *prob, double **Mu, double **LTSigma){
  int i, index_center, size_nb, *index_prob;
  int tmp_G = G - 1, tmp_min_n;
  double new_pi[1] = {1.0}, **new_Mu, **new_LTSigma, **new_X;
  double tmp_center;

  /* Get the seed state from R. */
  GetRNGstate();

  /* Use inverse CDF to sample a new center according to the given prob. */
  for(i = 1; i < n; i++) prob[i] = prob[i] + prob[i - 1];
  tmp_center = runif(0, prob[n - 1]); 

  if(tmp_center <= prob[0]){
    index_center = 0;
  } else{
    for(index_center = 1; index_center < n; index_center++){
      if(tmp_center > prob[index_center - 1] &&
         tmp_center <= prob[index_center]) break;
    }
  }

  /* Based on the new center to estimate the new ltsigma. */
  new_Mu = allocate_double_array(1);
  new_LTSigma = allocate_double_array(1);
  new_Mu[0] = Mu[tmp_G];
  new_LTSigma[0] = LTSigma[tmp_G];
  for(i = 0; i < p; i++) new_Mu[0][i] = X[index_center][i];
  est_ltsigma_mle_given_mu(X, n, p, new_Mu[0], new_LTSigma[0]);

  /* Compute prob based on the new center and ltsigma, and according
     to the prob to find the neighbors with size min.n + rpois(1, lambda). */
  for(i = 0; i < n; i++){
    prob[i] = mixllhd(p, 1, X[i], new_pi, new_Mu, new_LTSigma);
  }
  index_prob = (int *) orderDouble(prob, n); /* This is an increasing order. */
  size_nb = min_n + (int) rpois(lambda);

  /* Based on the neighbors to estimate Mu and LTSigma. */
  new_X = allocate_double_array(size_nb);
  tmp_min_n = n - size_nb;
  for(i = 0; i < size_nb; i++) new_X[i] = X[index_prob[tmp_min_n + i]];
  meandispersion_MLE(new_X, size_nb, p, new_Mu[0], new_LTSigma[0]);

  /* Release memory and set new seed state to R. */
  PutRNGstate();
  free(new_X);
  free(new_Mu);
  free(new_LTSigma);
  FREE_VECTOR(index_prob);
} /* End of cut_sub(). */


/* Model-based initializer. */
void mb_init(double **X, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma){
  int i, min_n = p + 1, n_nclass = (int) ceil((double) n / (double) nclass);
  int G, nc[nclass], i2, i3, j = 1;
  int tmp_G, tmp_class, class[n];
  double lambda, max_prob, prob[n], **gamm;
  double tmp_prob, tmp_n, tmp_pi, tmp_gamm;
  double class_prob[n], tmp_class_prob;
  PARAM param;
  void *pt_param;

  /* Find lambda, the expected size of neighbor. */
  if(n_nclass < min_n){
    //WCC printf("n is too small, or p or k is too large.\n"); 
    //WCC exit(1);
    error("n is too small, or p or k is too large.\n"); 
  }
  param.n_nclass = (double) (n_nclass - min_n);
  param.lower_bound = 1e-6;
  param.upper_bound = param.n_nclass;
  param.tol = 1e-6;
  param.maxit = 1000;
  pt_param = &param;
  lambda = find_lambda(pt_param);


  tmp_n = (double) n;
  MAKE_MATRIX(gamm, n, nclass);
  tmp_prob = 1.0 / (double) n;

  do{
    /* Initial prob. */
    for(i = 0; i < n; i++) prob[i] = tmp_prob;

    /* Note: cut_sub() will overwrite prob. */
    cut_sub(X, n, p, 1, min_n, lambda, prob, Mu, LTSigma);
    pi[0] = 1.0;
    for(i = 0; i < n; i++) class[i] = 0;
    nc[0] = n;

    for(G = 2; G <= nclass; G++){
      max_prob = 0.0;
      tmp_G = G - 1;

      for(i = 0; i < n; i++){
        prob[i] = mixllhd(p, tmp_G, X[i], pi, Mu, LTSigma);
        if(prob[i] > max_prob) max_prob = prob[i];
      }
      for(i = 0; i < n; i++) prob[i] = max_prob - prob[i];

      /* Drop the 75% points around cluster centers. */
      for(i = 0; i < G; i++){
        i3 = 0;
        for(i2 = 0; i2 < n; i2++){
          if(class[i2] == i) class_prob[i3++] = prob[i2];
        }
        tmp_class_prob = mb_quantile(i3, class_prob, 0.75);
        for(i2 = 0; i2 < n; i2++){
          if(class[i2] == i && prob[i2] < tmp_class_prob) prob[i2] = 0.0; 
        }
      }

      /* Note: cut_sub() will overwrite prob. */
      cut_sub(X, n, p, G, min_n, lambda, prob, Mu, LTSigma);

      /* Assume uniform for pi, do one estep, and reassign new pi. */
      tmp_pi = 1.0 / (double) G;
      for(i = 0; i < G; i++) pi[i] = tmp_pi;

      /* Run one estep to update pi. */
      estep(n, p, G, X, gamm, pi, Mu, LTSigma);
      for(i = 0; i < G; i++) nc[i] = 0;
      for(i = 0; i < n; i++){
        tmp_gamm = 0.0;
        tmp_class = 0;
        for(i2 = 0; i2 < G; i2++){
          if(tmp_gamm < gamm[i][i2]){
            tmp_gamm = gamm[i][i2];
            tmp_class = i2;
          }
        }
        class[i] = tmp_class;
        nc[tmp_class] = nc[tmp_class] + 1;
      }

      j = 1;
      for(i = 0; i < G; i++){
        if(nc[i] < min_n){
          j = 0;
	  break;
	}
      }
      if(j == 0) break;

      for(i = 0; i < G; i++) pi[i] = (double) nc[i] / tmp_n;
    }
  } while(j == 0);

  FREE_MATRIX(gamm);
} /* End of mb_init(). */

