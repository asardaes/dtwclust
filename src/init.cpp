#include <Rcpp.h>
#include "dtwclust.h"
#include "dtwclust++.h"

#define CALLDEF(name, n) { "C_"#name, (DL_FUNC) &name, n }
#define CALLDEFpp(name, n) { "C_"#name, (DL_FUNC) &dtwclust::name, n }

static R_CallMethodDef callMethods[] = {
    CALLDEFpp(SparseDistmatIndices__new, 1),
    CALLDEFpp(SparseDistmatIndices__getNewIndices, 4),
    CALLDEFpp(dba, 8),
    CALLDEFpp(dtwb_loop, 10),
    CALLDEFpp(dtw_lb, 5),
    CALLDEFpp(envelope, 2),
    CALLDEFpp(force_lb_symmetry, 1),
    CALLDEFpp(gak_loop, 9),
    CALLDEFpp(lbk, 4),
    CALLDEFpp(lbk_loop, 9),
    CALLDEFpp(lbi, 6),
    CALLDEFpp(lbi_loop, 11),
    CALLDEFpp(sbd_loop, 10),
    CALLDEFpp(tadpole, 8),
    CALLDEF(dtw_basic, 10),
    CALLDEF(logGAK, 8),
    CALLDEF(pairs, 2),
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
    DTWCLUST_REGISTER(pairs);
    DTWCLUST_REGISTER(setnames_inplace);
    DTWCLUST_REGISTER(tadpole);
    #undef DTWCLUST_REGISTER
}

extern "C" void R_init_dtwclust(DllInfo* info) {
    register_functions();
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
