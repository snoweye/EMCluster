#include <R.h>
#include <Rinternals.h>

SEXP R_init_EM(SEXP R_X, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_short_iter, SEXP R_short_eps, SEXP R_fixed_iter, SEXP R_n_candidate,
    SEXP R_EM_iter, SEXP R_EM_eps, SEXP R_lab, SEXP R_labK,
    SEXP R_init_method);
SEXP R_lnlikelihood(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_p_LTSigma, SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma);
SEXP R_mixllhd(SEXP R_x, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma);
SEXP R_dlmvnorm(SEXP R_x, SEXP R_p, SEXP R_p_LTsigma, SEXP R_mu,
    SEXP R_LTsigma);
SEXP R_emcluster(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_em_iter, SEXP R_em_eps);
SEXP R_estep(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_norm);
SEXP R_shortemcluster(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_p_LTSigma, SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma,
    SEXP R_maxiter, SEXP R_eps);
SEXP R_starts_via_svd(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_method, SEXP R_alpha);
SEXP R_assign(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma);
SEXP R_meandispersion(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_type);
SEXP R_M_emgroup(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_alpha, SEXP R_em_iter, SEXP R_em_eps);
SEXP R_mstep(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_Gamma);
SEXP R_RRand(SEXP R_N, SEXP R_TRUK, SEXP R_PREDK, SEXP R_trcl, SEXP R_prcl);
SEXP ss_R_emcluster(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass,
    SEXP R_p_LTSigma, SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_em_iter,
    SEXP R_em_eps, SEXP R_lab);
SEXP ss_R_assign(SEXP R_x, SEXP R_n, SEXP R_p, SEXP R_nclass, SEXP R_p_LTSigma,
    SEXP R_pi, SEXP R_Mu, SEXP R_LTSigma, SEXP R_lab);
