/* This file contains function for initializing Mus and LTSigmas.
   Modified: Wei-Chen Chen 2009/03/14.
*/

#include "array.h"

/* This function is provided by Volodymyr Melnykov. */
void labInitMus(int n, int p, int labK, double **x, int *lab, double **Mu){
  int i, j;
  double *labNk;

  MAKE_VECTOR(labNk, labK);
  
  for(i = 0; i < labK; i++){
    labNk[i] = 0;
    for(j = 0; j < p; j++) Mu[i][j] = 0;
  }
  
  for(i = 0; i < n; i++){
    if(lab[i] != -1){
      for(j = 0; j < p; j++) Mu[lab[i]][j] = Mu[lab[i]][j] + x[i][j];
      labNk[lab[i]] = labNk[lab[i]] + 1;
    }
  }

  for(i = 0; i < labK; i++){
    for(j = 0; j < p; j++) Mu[i][j] = Mu[i][j] / labNk[i];
  }

  FREE_VECTOR(labNk);
} /* End of labInitMUs(). */
