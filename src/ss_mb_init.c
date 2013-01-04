/* This file is modified from "mb_init.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"


/* Model-based initializer. */
void ss_mb_init(double **X, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma, int *lab, int labK){
  int i, min_n = p + 1, n_nclass = (int) ceil((double) n / (double) nclass);
  int G, nc[nclass], i2, i3, j = 1;
  int tmp_G, tmp_class = 0, class[n];
  double lambda, max_prob, prob[n], **gamm;
  double tmp_n, tmp_pi, tmp_gamm;
  double class_prob[n], tmp_class_prob;
  PARAM param;
  void *pt_param;

  int lab_n = 0;
  double **lab_X, **lab_gamm;

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
  for(i = 0; i < n; i++){
    if(lab[i] > -1) lab_n = lab_n + 1;
  }
  lab_X = allocate_double_array(lab_n);
  lab_gamm = allocate_double_array(lab_n);
  i3 = 0;
  for(i = 0; i < n; i++){
    if(lab[i] > -1){
      lab_X[i3] = X[i];
      lab_gamm[i3++] = gamm[i];
    }
  }

  do{
    /* Initial for known clusters. */
    for(i = 0; i < n; i++){
      if(lab[i] > -1){
        for(i2 = 0; i2 < nclass; i2++) gamm[i][i2] = 0.0;
        gamm[i][lab[i]] = 1.0;
      }
    }
    mstep(lab_X, lab_n, p, labK, pi, Mu, LTSigma, lab_gamm);

    for(G = labK + 1; G <= nclass; G++){
      max_prob = 0.0;
      tmp_G = G - 1;

      for(i = 0; i < n; i++){
        prob[i] = mixllhd(p, tmp_G, X[i], pi, Mu, LTSigma);
        if(prob[i] > max_prob) max_prob = prob[i];
      }
      for(i = 0; i < n; i++) prob[i] = max_prob - prob[i];

      /* Drop the 75% points around cluster centers. */
      for(i = labK; i < G; i++){
        i3 = 0;
        for(i2 = 0; i2 < n; i2++){
          if(class[i2] == i) class_prob[i3++] = prob[i2];
        }
        tmp_class_prob = mb_quantile(i3, class_prob, 0.75);
        for(i2 = 0; i2 < n; i2++){
          if(class[i2] == i && prob[i2] < tmp_class_prob) prob[i2] = 0.0; 
        }
      }

      /* Set prob = 0 for known clusters where to avoid pick a center from. */
      for(i = 0; i < n; i++){
        if(lab[i] == -1) prob[i] = 0.0;
      }

      /* Note: ss_cut_sub() will overwrite prob. */
      cut_sub(X, n, p, G, min_n, lambda, prob, Mu, LTSigma);

      /* Assume uniform for pi, do one estep, and reassign new pi. */
      tmp_pi = 1.0 / (double) G;
      for(i = 0; i < G; i++) pi[i] = tmp_pi;

      /* Run one estep to update pi. */
      ss_estep(n, p, G, X, gamm, pi, Mu, LTSigma, lab);
      for(i = 0; i < G; i++) nc[i] = 0;
      for(i = 0; i < n; i++){
        if(lab[i] == -1){
          tmp_gamm = 0.0;
          tmp_class = 0;
          for(i2 = 0; i2 < G; i2++){
            if(tmp_gamm < gamm[i][i2]){
              tmp_gamm = gamm[i][i2];
              tmp_class = i2;
            }
          }
          class[i] = tmp_class;
        } else{
          class[i] = lab[i];
        }
        nc[tmp_class] = nc[tmp_class] + 1;
      }

      j = 1;
      for(i = 0; i < G; i++){
        if(nc[i] <= p){
          j = 0;
	  break;
	}
      }
      if(j == 0) break;

      for(i = 0; i < G; i++) pi[i] = (double) nc[i] / tmp_n;
    }
  } while(j == 0);

  FREE_MATRIX(gamm);
  free(lab_X);
  free(lab_gamm);
} /* End of ss_mb_init(). */

