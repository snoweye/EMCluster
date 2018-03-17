/* This file is modified from "init_other.c" for semi-supervised clustering.
   Required int *lab,
     -1 for clusters unknown,
     0 to (labK - 1) for data with known "labK" clusters.
   Modified: Wei-Chen Chen 2009/03/12.
*/

#include "ss_mb_tool.h"
#include "mat_vec.h"

/* This function is called by ss_shortems().
   Mu[0, ..., labK-1] should be assigned before calling this function.
*/
void ss_randomEMinit(double **x, int n, int p, int nclass, double *pi,
   double **Mu, double **LTSigma,
   int *lab, int labK, int nonlab_total, int *lab_index){
  int *ordr, i, j, *clas, *nc;
  int new_nclass = nclass - labK;
  double labMu[labK][p];
  
  for(i = 0; i < labK; i++){
    for(j = 0; j < p; j++) labMu[i][j] = Mu[i][j];
  }
  
  /* Initial centers for all other unknown clusters. */
  MAKE_VECTOR(ordr, new_nclass);
  MAKE_VECTOR(clas, n);
  MAKE_VECTOR(nc, nclass);
  do{
    for(i = 0; i < labK; i++){
      for(j = 0; j < p; j++) Mu[i][j] = labMu[i][j];
    }
    i = srswor(nonlab_total, new_nclass, ordr);
    for(i = labK; i < nclass; i++){
      for(j = 0; j < p; j++) Mu[i][j] = x[lab_index[ordr[i - labK]]][j];
    }
    for(i = 0; i < n; i++){
      if(lab[i] == -1){
        clas[i] = assign_closest(x[i], p, nclass, Mu);
      } else{
        clas[i] = lab[i];
      }
    } 
    j = initials(x, n, p, nclass, nc, Mu, LTSigma, clas);
  } while(j == 0);
  for(i = 0; i < nclass; i++) pi[i] = 1. * nc[i] / n;
  FREE_VECTOR(nc);
  FREE_VECTOR(clas);
  FREE_VECTOR(ordr);
} /* End of ss_randomEMinit(). */


int ss_shortemcluster_org(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int maxiter, double eps, double *llhdval,
    int *lab){
  int iter;
  double **gamm, llhd, oldllhd, llh0;
  /*same as emcluster, only difference being in how convergence is handled,
   done as per Biernacki et al*/
  MAKE_MATRIX(gamm, n, k);
  llhd = lnlikelihood(n, p, k, pi, X, Mu, LTSigma);
  llh0 = llhd;
  iter = 0;
  do{
    oldllhd = llhd;
    ss_estep(n, p, k, X, gamm, pi, Mu, LTSigma, lab);
    mstep(X, n, p, k, pi, Mu, LTSigma, gamm);
    llhd = lnlikelihood(n, p, k, pi, X, Mu, LTSigma);
    iter++;
  } while(((oldllhd - llhd) / (llh0 - llhd) > eps) && (iter < maxiter));
/* Marked by Wei-Chen Chen on 2009/01/21.
   This value, (oldllhd-llhd)/(llhd-llh0), is always negative.

  } while ((((oldllhd-llhd)/(llhd-llh0)) > eps) && (iter<maxiter));
*/
  *llhdval = llhd;
  FREE_MATRIX(gamm);
  return iter;
} /* End of ss_shortemcluster_org(). */


/* Created by Wei-Chen Chen on 2009/06/16. */
int ss_shortemcluster(int n, int p, int k, double *pi, double **X,
    double **Mu, double **LTSigma, int maxiter, double eps, double *llhdval,
    int *conv_iter, double *conv_eps, int *lab){
  int iter, i, n_par =  p * (p + 1) / 2;
  double *backup_pi, **backup_Mu, **backup_LTSigma;
  double **gamm, llhd, oldllhd, llh0;

  MAKE_VECTOR(backup_pi, k);
  MAKE_MATRIX(backup_Mu, k, p);
  MAKE_MATRIX(backup_LTSigma, k, n_par);
  MAKE_MATRIX(gamm, n, k);

  estep_gamma(n, p, k, X, gamm, Mu, LTSigma);
  llhd = lnlikelihood_gamma(n, k, gamm, pi);
  llh0 = llhd;
  iter = 0;
  do{
    oldllhd = llhd;
    ss_norm_gamma(n, k, gamm, pi, lab);

    for(i = 0; i < k; i++) backup_pi[i] = pi[i];
    cpy(Mu, k, p, backup_Mu);
    cpy(LTSigma, k, n_par, backup_LTSigma);
    
    mstep(X, n, p, k, pi, Mu, LTSigma, gamm);
    estep_gamma(n, p, k, X, gamm, Mu, LTSigma);
    llhd = lnlikelihood_gamma(n, k, gamm, pi);

    if(oldllhd > llhd){
      for(i = 0; i < k; i++) pi[i] = backup_pi[i];
      cpy(backup_Mu, k, p, Mu);
      cpy(backup_LTSigma, k, n_par, LTSigma);
      llhd = oldllhd;
      iter--;
      break;
    }
    
    iter++;
    *conv_eps = fabs((oldllhd - llhd) / (llh0 - llhd));
  } while((*conv_eps > eps) && (iter < maxiter));
  *llhdval = llhd;
  *conv_iter = iter;

  FREE_VECTOR(backup_pi);
  FREE_MATRIX(backup_Mu);
  FREE_MATRIX(backup_LTSigma);
  FREE_MATRIX(gamm);
  return iter;
} /* End of ss_shortemcluster(). */


int ss_shortems(int n, int p, int nclass, double *pi, double **X, double **Mu,  
    double **LTSigma, int maxshortiter, double shorteps,
    int *conv_iter, double *conv_eps, int *lab, int labK){
  /*initializing as per Beiernacki, Celeaux, Govaert~(2003) */

  int i, j, iter, totiter = 0, n_par = p * (p + 1) / 2;
  int nonlab_total = 0, lab_index[n];
  double *oldpi, **oldMu, **oldLTSigma, oldllh = -Inf, llhval;
  double **labMu;

  MAKE_VECTOR(oldpi, nclass);
  MAKE_MATRIX(oldMu, nclass, p);
  MAKE_MATRIX(oldLTSigma, nclass, n_par);
  MAKE_MATRIX(labMu, labK, p);

  for(i = 0; i < n; i++){
    if(lab[i] == -1) lab_index[nonlab_total++] = i;
  }
  labInitMus(n, p, labK, X, lab, labMu);

  do{
    for(i = 0; i < labK; i++){
      for(j = 0; j < p; j++) oldMu[i][j] = labMu[i][j];
    }

    iter = maxshortiter - totiter;

/* Modified by Wei-Chen Chen on 2009/03/08.
    ss_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma,
                    lab, labK, nonlab_total, lab_index);
    ss_mb_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma,
                       lab, labK, nonlab_total, lab_index);
*/
    ss_randomEMinit(X, n, p, nclass, oldpi, oldMu, oldLTSigma,
                    lab, labK, nonlab_total, lab_index);

    iter = ss_shortemcluster(n, p, nclass, oldpi, X, oldMu, oldLTSigma, iter,
                             shorteps, &llhval, conv_iter, conv_eps, lab);
    if(llhval >= oldllh){
      int i;
      oldllh = llhval;
      cpy(oldMu, nclass, p, Mu);
      cpy(oldLTSigma, nclass, n_par, LTSigma);
      for(i = 0; i < nclass; i++) pi[i] = oldpi[i];
    }
    totiter += iter;
  } while(totiter < maxshortiter);
  FREE_MATRIX(oldMu);
  FREE_MATRIX(oldLTSigma);
  FREE_VECTOR(oldpi);
  FREE_MATRIX(labMu);
  return totiter; 
} /* End of ss_shortems(). */


/*
  This routine clusters a dataset into an optimal number of groups using the 
  E-M algorithm and a choice of AIC or BIC and with start points provided by 
  the shortems function.  
  The function returns    the classification ids as well as the parameter 
  estimates. The parameter estimates are created in the routine and should not
  be allocated before the function call. This also means that the parameter 
  estimates have to be cleaned after the function call. (NOTE: This is an 
  important fact to remember. */
int ss_em_EM(double **x, int n, int p, int nclass, double *pi, double **Mu,
    double **LTSigma, double *llhdval, int *nc, int shortiter,
    double shorteps, int *conv_iter, double *conv_eps, int *lab, int labK){
  int flag = 0;
  double like;

  if(nclass == 1){
/* These formulae are not correct with two errors,
   1. n-1 should be replace by n in meandispersion() for LTSigma, and
   2. determinant() will change the values of LTSigma[0], see "initials.c".
   Modified: Wei-Chen Chen on 2008/12/05.

    meandispersion(x,n,p,Mu[0],LTSigma[0]);
*/
    meandispersion_MLE(x, n, p, Mu[0], LTSigma[0]);
    like = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
           0.5 * n * p * log(2 * PI);
    *llhdval = like;
  }
  else {
    ss_shortems(n, p, nclass, pi, x, Mu, LTSigma, shortiter, shorteps,
                conv_iter, conv_eps, lab, labK);
    ss_emcluster(n, p, nclass, pi, x, Mu, LTSigma, 1000, 0.0001, &like,
                 conv_iter, conv_eps, lab);
    *llhdval = like;
  } 
  return flag;
} /* End of ss_em_EM(). */

