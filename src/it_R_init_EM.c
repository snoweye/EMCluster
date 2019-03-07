/* This file contains functions called by R wraps in "R/fcn_init_EM.r",
   and these functions call the relative functions in "src/it_*.c"
   using in Dr. Maitra's ten-clusters programs.

   Writen: Wei-Chen Chen on 2009/04/29.

   Move PROTECT from the create_emptr().
   Modified: Wei-Chen Chen on 2019/03/07.
*/

#include "it_tool.h"

/* This function calls other initializations and is called by
   init.EM() using .Call() in "R/fcn_init_EM.r".
   Input:
     R_X: SEXP[R_n * R_p], data matrix of R_n*R_p.
     R_n: SEXP[1], number of observations.
     R_p: SEXP[1], number of dimersions.
     R_nclass: SEXP[1], number of classes.
     R_short_iter: SEXP[1], number of short iterations, 500 by default.
     R_short_eps: SEXP[1], epsilon for short iterations, 1e-2 by default.
     R_fixed_iter: SEXP[1], number of rand iterations, 1 by default.
     R_n_candidate: SEXP[1], number of candidate, 3 by default.
     R_EM_iter: SEXP[1], max iterations for emclust(), 1000 by default.
     R_EM_eps: SEXP[1], tolerance for emclust(), 1e-4 by default.
     R_lab: SEXP[1], -1 for points with unknown clusters;
                     0, .., (labK-1) for known.
     R_labK: SEXP[1], the number of known clusters.
     R_init_method: SEXP[1], initialization method. (see below for detail)
#####
#   unsupervised init_method:
#     1 = em.EM,     2 = Rnd.EM,
#     11 = MBem.EM, 12 = MBRnd.EM
#     21 = acem.EM, 22 = acRnd.EM
#   semi-supervised init_method:
#     101 = em.EM,   102 = Rnd.EM
#     111 = MBem.EM, 112 = MBRnd.EM
#   For iteration & time record:
#     3 = oneRnd,      13 = oneMBRnd
#     103 = ss.oneRnd, 113 = ss.oneMBRnd
#####
   Output:
     emobj: a list contains
       pi: SEXP[R_nclass], proportions of classes.
       Mu: SEXP[R_nclass, R_p], means of MVNs.
       LTSigma: SEXP[nclass, R_p * (R_p + 1) / 2], lower triangular
                sigma matrices.
       llhdval: SEXP[1], log likelihood value.
       nc: SEXP[R_nclass], number of observations in each class.
       class: SEXP[R_n], class id's for all observations
              starting from 0 to (R_nclass - 1).
       conv_iter: SEXP[1], convergent iterations.
       conv_eps: SEXP[1], convergent tolerance.
*/
EMPTR allocate_emptr(void){
  EMPTR emptr = (EMPTR) malloc(sizeof(struct emptr));
  if(emptr == NULL){
    error("Memory allocation fails!\n");
  }
  return(emptr);
}

SEXP R_init_EM(SEXP R_X, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_short_iter, SEXP R_short_eps, SEXP R_fixed_iter, SEXP R_n_candidate,
    SEXP R_EM_iter, SEXP R_EM_eps, SEXP R_lab, SEXP R_labK,
    SEXP R_init_method){
  /* int C_protect_length; */
  EMPTR emptr = allocate_emptr();
  SEXP emobj;

  /* Initial emptr. */
  PROTECT(emobj = create_emptr(R_X, R_n, R_p, R_nclass,
                               R_short_iter, R_short_eps, R_fixed_iter,
                               R_n_candidate,
                               R_EM_iter, R_EM_eps, R_lab, R_labK,
                               R_init_method, emptr));

  /* Compute. */
  init_EM(emptr->C_X, emptr->C_pi, emptr->C_Mu, emptr->C_LTSigma,
          emptr->C_llhdval,
          *emptr->C_n, *emptr->C_p, *emptr->C_nclass, emptr->C_nc,
          emptr->C_class,
          *emptr->C_short_iter, *emptr->C_short_eps, *emptr->C_fixed_iter,
	  *emptr->C_n_candidate,
          *emptr->C_EM_iter, *emptr->C_EM_eps,
          emptr->C_conv_iter, emptr->C_conv_eps,
          emptr->C_lab, *emptr->C_labK,
          *emptr->C_init_method);

  /* Free memory and release protectation. */
  free(emptr->C_X);
  free(emptr->C_Mu);
  free(emptr->C_LTSigma);
  /* C_protect_length = emptr->C_protect_length; */
  free(emptr);

  /* UNPROTECT(C_protect_length); */
  UNPROTECT(1);
  return(emobj);
} /* End of R_init_EM(). */

