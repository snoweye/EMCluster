/* This file contains several functions to estimate mean and dispersion
   including MLE and MME methods modified from meandispersion() in
   "initials.c" provided by Dr. Maitra.
   Created: Wei-Chen Chen 2009/02/11.
*/


/* This function provides a correct ltsigma for MLE.
   Modified: Wei-Chen Chen on 2008/12/05.
*/
void meandispersion_MLE(double **x, int n, int p, double *mu, double *ltsigma){
  /* This routine calculates the mean and dispersion of a homogeneous sample */
  int i, j, l, n_par = p * (p + 1) / 2, tmp_int[p];
  double tmp_double;

  for(i = 0; i < p; i++) mu[i] = 0.;
  for(i = 0; i < n; i++){
    for(j = 0; j < p; j++) mu[j] += x[i][j];
  }
  for(j = 0; j < p; j++) mu[j] /= n;

  for(i = 0; i < n_par; i++) ltsigma[i] = 0.;
  for(j = 0; j < p; j++) tmp_int[j] = j * (j + 1) / 2;
  for(i = 0; i < n; i++){
    for(j = 0; j < p; j++){
      for(l = 0; l <= j; l++)
        ltsigma[tmp_int[j] + l] += (x[i][j] - mu[j]) * (x[i][l] - mu[l]);
    }
  }
  if(n > 1) {
    tmp_double = (double) n;
    for(j = 0; j < n_par; j++) ltsigma[j] /= tmp_double;
  }
} /* End of meandispersion_MLE(). */


/* This function provides a correct ltsigma for MME.
   Modified: Wei-Chen Chen on 2008/12/05.
*/
void meandispersion_MME(double **x, int n, int p, double *mu, double *ltsigma){
  /* This routine calculates the mean and dispersion of a homogeneous sample */
  int i,j,l, n_par = p * (p + 1) / 2, tmp_int[p];
  double tmp_double;

  for(i = 0; i < p; i++) mu[i] = 0.;
  for(i = 0; i < n; i++){
    for(j = 0; j < p; j++) mu[j] += x[i][j];
  }
  for(j = 0; j < p; j++) mu[j] /= n;

  for(i = 0; i < n_par; i++) ltsigma[i] = 0.0;
  for(j = 0; j < p; j++) tmp_int[j] = j * (j + 1) / 2;
  for(i = 0; i < n; i++){
    for(j = 0; j < p; j++){
      for(l = 0; l <= j; l++)
        ltsigma[tmp_int[j] + l] += (x[i][j] - mu[j]) * (x[i][l] - mu[l]);
    }
  }
  if(n > 1){
    tmp_double = (double) (n - p);
    for(j = 0; j < n_par; j++) ltsigma[j] /= tmp_double;
  }
}


/* This function provides a ltsigma for MLE by given mu.
   Modified: Wei-Chen Chen on 2009/02/10.
*/
void est_ltsigma_mle_given_mu(double **x, int n, int p, double *mu,
    double *ltsigma){
  /* This routine calculates the mean and dispersion of a homogeneous sample */
  int i, j, l, n_par = p * (p + 1) / 2, tmp_int[p];
  double tmp_double;

  for(i = 0; i < n_par; i++) ltsigma[i] = 0.0;
  for(j = 0; j < p; j++) tmp_int[j] = j * (j + 1) / 2;
  for(i = 0; i < n; i++){
    for(j = 0; j < p; j++){
      for(l = 0; l <= j; l++){
        ltsigma[tmp_int[j] + l] += (x[i][j] - mu[j]) * (x[i][l] - mu[l]);
      }
    }
  }
  if(n > 1){
    tmp_double = (double) n;
    for(j = 0; j < n_par; j++) ltsigma[j] /= tmp_double;
  }
}
