/* This file contains functions called by R wraps in "src/it_R_init_EM.r".
   using in Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2009/04/29.
*/

#include "it_tool.h"

void init_EM(double **X, double *pi, double **Mu, double **LTSigma,
    double *llhdval,
    int n, int p, int nclass, int *nc, int *class,
    int short_iter, double short_eps, int fixed_iter, int n_candidate,
    int EM_iter, double EM_eps, int *conv_iter, double *conv_eps,
    int *lab, int labK,
    int init_method){
  int j;
//  clock_t start;

  if(nclass == 1){
    nc[0] = n;
    pi[0] = 1.0;
    for(j = 0; j < n; j++) class[j] = 0;
    meandispersion_MLE(X, n, p, Mu[0], LTSigma[0]);
    *llhdval = -0.5 * n * p - 0.5 * n * log(determinant(LTSigma[0], p)) -
               0.5 * n * p * log(2 * PI);
  } else {
/* Unsupervised clustering. */
    if(init_method == 1){	/* em.EM */
      shortems(n, p, nclass, pi, X, Mu, LTSigma,
               short_iter, short_eps);
    } else if(init_method == 2){	/* Rnd.EM */
      mod_shortems(n, p, nclass, pi, X, Mu, LTSigma,
                   short_iter, fixed_iter);
/* Model-based unsupervised clustering. */
    } else if(init_method == 11){	/* MBem.EM */
      shortems_mb(n, p, nclass, pi, X, Mu, LTSigma,
                  short_iter, short_eps);
    } else if(init_method == 12){	/* MBRnd.EM */
      mod_shortems_mb(n, p, nclass, pi, X, Mu, LTSigma,
                      short_iter, fixed_iter);
/* Adapted candidated unsupervised clustering. */
//    } else if(init_method == 21){	/* acem.EM */
//      shortems_ac(n, p, nclass, pi, X, Mu, LTSigma,
//                  short_iter, short_eps, n_candidate);
//    } else if(init_method == 22){	/* acRnd.EM */
//      mod_shortems_ac(n, p, nclass, pi, X, Mu, LTSigma,
//                      short_iter, fixed_iter, n_candidate);
/* Semi-supervised clustering. */
    } else if(init_method == 101){	/* ss.em.EM */
      ss_shortems(n, p, nclass, pi, X, Mu, LTSigma,
                  short_iter, short_eps,
                  lab, labK);
    } else if(init_method == 102){	/* ss.Rnd.EM */
      ss_mod_shortems(n, p, nclass, pi, X, Mu, LTSigma,
                      short_iter, fixed_iter,
                      lab, labK);
/* Model-based semi-supervised clustering. */
    } else if(init_method == 111){	/* ss.MBem.EM */
      ss_shortems_mb(n, p, nclass, pi, X, Mu, LTSigma,
                     short_iter, short_eps,
                     lab, labK);
    } else if(init_method == 112){	/* ss.MBRnd.EM */
      ss_mod_shortems_mb(n, p, nclass, pi, X, Mu, LTSigma,
                         short_iter, fixed_iter,
                         lab, labK);
/* Adapted candidated semi-supervised clustering. */
//    } else if(init_method == 121){	/* ss.em.EM */
//      ss_shortems_ac(n, p, nclass, pi, X, Mu, LTSigma,
//                     short_iter, short_eps, n_candidate,
//                     lab, labK);
//    } else if(init_method == 122){	/* ss.Rnd.EM */
//      ss_mod_shortems_ac(n, p, nclass, pi, X, Mu, LTSigma,
//                         short_iter, fixed_iter, n_candidate,
//                         lab, labK);
/* For output detail iterations. */
//    } else if(init_method == -1){	/* oneRnd */
//      Rprintf("iter\tdt\tlogL\n");
//      start = clock();
//      randomEMinit(X, n, p, nclass, pi, Mu, LTSigma);
//      Rprintf("%d\t%16.10lf\t%24.16lf\n", -1,
//              (double) (clock() - start) / (double) CLOCKS_PER_SEC, 0);
//    } else if(init_method == -11){	/* oneMBRnd */
//      Rprintf("iter\tdt\tlogL\n");
//      start = clock();
//      mb_init(X, n, p, nclass, pi, Mu, LTSigma);
//      Rprintf("%d\t%16.10lf\t%24.16lf\n", -1,
//              (double) (clock() - start) / (double) CLOCKS_PER_SEC, 0);
//    } else if(init_method == -101){	/* ss.oneRnd */
//      Rprintf("iter\tdt\tlogL\n");
//      start = clock();
//
//      int i, j, nonlab_total = 0, lab_index[n];
//      double **labMu;
//
//      MAKE_MATRIX(labMu, labK, p);
//      for(i = 0; i < n; i++){
//        if(lab[i] == -1) lab_index[nonlab_total++] = i;
//      }
//      labInitMus(n, p, labK, X, lab, labMu);
//      for(i = 0; i < labK; i++){
//        for(j = 0; j < p; j++) Mu[i][j] = labMu[i][j];
//      }
//      free(labMu);
//
//      ss_randomEMinit(X, n, p, nclass, pi, Mu, LTSigma,
//                         lab, labK, nonlab_total, lab_index);
//      Rprintf("%d\t%16.10lf\t%24.16lf\n", -1,
//              (double) (clock() - start) / (double) CLOCKS_PER_SEC, 0);
//    } else if(init_method == -111){	/* ss.oneMBRnd */
//      Rprintf("iter\tdt\tlogL\n");
//      start = clock();
//      ss_mb_init(X, n, p, nclass, pi, Mu, LTSigma, lab, labK);
//      Rprintf("%d\t%16.10lf\t%24.16lf\n", -1,
//              (double) (clock() - start) / (double) CLOCKS_PER_SEC, 0);
    } else{
      error("Method is not found.");
    }

/* Unsupervised clustering. */
    if(init_method == 1 || init_method == 2 ||
       init_method == 11 || init_method == 12){
//       init_method == 11 || init_method == 12 ||
//       init_method == 21 || init_method == 22){
      emcluster(n, p, nclass, pi, X, Mu, LTSigma, EM_iter, EM_eps, llhdval,
                conv_iter, conv_eps);
      assign(n, p, nclass, X, pi, Mu, LTSigma, class, nc);
/* Semi-supervised clustering. */
    } else if(init_method == 101 || init_method == 102 ||
              init_method == 111 || init_method == 112){
//              init_method == 111 || init_method == 112 ||
//              init_method == 121 || init_method == 122){
      ss_emcluster(n, p, nclass, pi, X, Mu, LTSigma, EM_iter, EM_eps, llhdval,
                   lab);
      ss_assign(n, p, nclass, X, pi, Mu, LTSigma, class, nc, lab);
/* For output detail iterations. */
//    } else if(init_method == -1 || init_method == -11){
//      it_emcluster(n, p, nclass, pi, X, Mu, LTSigma, EM_iter, EM_eps, llhdval);
//      assign(n, p, nclass, X, pi, Mu, LTSigma, class, nc);
//    } else if(init_method == -101 || init_method == -111){
//      it_ss_emcluster(n, p, nclass, pi, X, Mu, LTSigma, EM_iter, EM_eps,
//                      llhdval, lab);
//      ss_assign(n, p, nclass, X, pi, Mu, LTSigma, class, nc, lab);
    } else{
      error("Method is not found.");
    }
  } 

} /* END of init_EM(). */

