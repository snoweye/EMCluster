/* This file is modified from "initials.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/13.
*/

#include "ss_mb_tool.h"

void ss_assign(int n, int p, int k, double **X, double *pi, double **Mu,
    double **LTSigma, int *class, int *nc, int *lab){
  int i;

  for(i = 0; i < k; i++)  nc[i] = 0;
  for(i = 0; i < n; i++){
    if(lab[i] == -1){
      class[i] = classify(X[i], p, k, pi, Mu, LTSigma);
    } else{
      class[i] = lab[i];
    }
    nc[class[i]]++;
  }
}
