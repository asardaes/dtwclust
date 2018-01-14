#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

#include "dtwclust++.h"
#include "dtwclust.h"

#define CALLDEF(name, n) { "C_"#name, (DL_FUNC) &name, n }
#define CALLDEFpp(name, n) { "C_"#name, (DL_FUNC) &dtwclust::name, n }

static R_CallMethodDef callMethods[] = {
    CALLDEFpp(SparseDistmatIndices__new, 1),
    CALLDEFpp(SparseDistmatIndices__getNewIndices, 4),
    CALLDEFpp(dba, 9),
    CALLDEFpp(distmat_loop, 8),
    CALLDEFpp(dtw_lb, 6),
    CALLDEFpp(envelope, 2),
    CALLDEFpp(force_lb_symmetry, 1),
    CALLDEFpp(lbk, 4),
    CALLDEFpp(lbi, 6),
    CALLDEFpp(sdtw_cent, 8),
    CALLDEFpp(soft_dtw, 6),
    CALLDEFpp(tadpole, 9),
    CALLDEF(dtw_basic, 11),
    CALLDEF(logGAK, 8),
    CALLDEF(pairs, 1),
    CALLDEF(setnames_inplace, 2),
    {NULL, NULL, 0}
};

void register_functions() {
    using namespace dtwclust;

    #define DTWCLUST_REGISTER(__FUN__) R_RegisterCCallable("dtwclust", #__FUN__, (DL_FUNC)__FUN__)
    DTWCLUST_REGISTER(dba);
    DTWCLUST_REGISTER(dtw_basic);
    DTWCLUST_REGISTER(dtw_lb);
    DTWCLUST_REGISTER(envelope);
    DTWCLUST_REGISTER(lbk);
    DTWCLUST_REGISTER(lbi);
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
