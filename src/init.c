#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP pairs(SEXP, SEXP);
SEXP envelop(SEXP, SEXP);
SEXP dtw_basic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP logGAK(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
    {"C_pairs", (DL_FUNC) &pairs, 2},
    {"C_envelop", (DL_FUNC) &envelop, 2},
    {"C_dtw_basic", (DL_FUNC) &dtw_basic, 10},
    {"C_logGAK", (DL_FUNC) &logGAK, 8},
    {NULL, NULL, 0}
};


void R_init_dtwclust(DllInfo* info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
