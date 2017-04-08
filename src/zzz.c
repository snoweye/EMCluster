#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "it_tool.h"

static const R_CallMethodDef callMethods[] = {
	{"create_emptr", (DL_FUNC) &create_emptr, 14},
	{"R_init_EM", (DL_FUNC) &R_init_EM, 13},
	{"R_lnlikelihood", (DL_FUNC) &R_lnlikelihood, 8},
	{"R_mixllhd", (DL_FUNC) &R_mixllhd, 7},
	{"R_dlmvnorm", (DL_FUNC) &R_dlmvnorm, 5},
	{"R_lnlikelihood", (DL_FUNC) &R_lnlikelihood, 8},
	{"R_mixllhd", (DL_FUNC) &R_mixllhd, 7},
	{"R_dlmvnorm", (DL_FUNC) &R_dlmvnorm, 5},
	{"R_emcluster", (DL_FUNC) &R_emcluster, 10},
	{"R_estep", (DL_FUNC) &R_estep, 9},
	{"R_shortemcluster", (DL_FUNC) &R_shortemcluster, 10},
	{"R_starts_via_svd", (DL_FUNC) &R_starts_via_svd, 6},
	{"R_assign", (DL_FUNC) &R_assign, 8},
	{"R_meandispersion", (DL_FUNC) &R_meandispersion, 4},
	{"R_M_emgroup", (DL_FUNC) &R_M_emgroup, 7},
	{"R_mstep", (DL_FUNC) &R_mstep, 5},
	{"R_RRand", (DL_FUNC) &R_RRand, 5},
	{"ss_R_emcluster", (DL_FUNC) &ss_R_emcluster, 11},
	{"ss_R_assign", (DL_FUNC) &ss_R_assign, 9},

	/* Finish R_CallMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_EMCluster(DllInfo *info){
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_EMCluster(). */
