#include <Rcpp.h>
#include "dtwclust.h"
#include "dtwclustpp.h"

#define CALLDEF(name, n) { "C_"#name, (DL_FUNC) &name, n }
#define CALLDEFpp(name, n) { "C_"#name, (DL_FUNC) &dtwclust::name, n }

static R_CallMethodDef callMethods[] = {
    CALLDEFpp(dba, 8),
    CALLDEFpp(envelope, 2),
    CALLDEFpp(lbk, 4),
    CALLDEFpp(lbi, 6),
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
    DTWCLUST_REGISTER(envelope);
    DTWCLUST_REGISTER(lbk);
    DTWCLUST_REGISTER(lbi);
    DTWCLUST_REGISTER(logGAK);
    DTWCLUST_REGISTER(pairs);
    DTWCLUST_REGISTER(setnames_inplace);
    #undef DTWCLUST_REGISTER
}

extern "C" void R_init_dtwclust(DllInfo* info) {
    register_functions();
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
