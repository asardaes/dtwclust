#include "dtwclust.h"

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n) { "C_"#name, (DL_FUNC) &dtwclust::name, n }

static R_CallMethodDef callMethods[] = {
    CALLDEF(SparseDistmatIndices__new, 1),
    CALLDEF(SparseDistmatIndices__getNewIndices, 4),
    CALLDEF(dba, 9),
    CALLDEF(distmat_loop, 8),
    CALLDEF(dtw_basic, 11),
    CALLDEF(dtw_lb, 6),
    CALLDEF(envelope, 2),
    CALLDEF(force_lb_symmetry, 1),
    CALLDEF(lbk, 4),
    CALLDEF(lbi, 6),
    CALLDEF(logGAK, 8),
    CALLDEF(pairs, 1),
    CALLDEF(setnames_inplace, 2),
    CALLDEF(sdtw_cent, 6),
    CALLDEF(soft_dtw, 5),
    CALLDEF(tadpole, 9),
    {NULL, NULL, 0}
};

void register_functions() {
    using namespace dtwclust;

    #define DTWCLUST_REGISTER(__FUN__) R_RegisterCCallable("dtwclust", #__FUN__, (DL_FUNC)__FUN__)
    DTWCLUST_REGISTER(dba);
    DTWCLUST_REGISTER(dtw_basic);
    DTWCLUST_REGISTER(dtw_lb);
    DTWCLUST_REGISTER(envelope);
    DTWCLUST_REGISTER(lbi);
    DTWCLUST_REGISTER(lbk);
    DTWCLUST_REGISTER(logGAK);
    DTWCLUST_REGISTER(sdtw_cent);
    DTWCLUST_REGISTER(soft_dtw);
    DTWCLUST_REGISTER(tadpole);
    #undef DTWCLUST_REGISTER
}

extern "C" void R_init_dtwclust(DllInfo* info) {
    register_functions();
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
