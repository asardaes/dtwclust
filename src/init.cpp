#include <Rcpp.h>
#include "dtwclust.h"
#include "dtwclustpp.h"

#define CALLDEF(name, n) { "C_"#name, (DL_FUNC) &name, n }

static R_CallMethodDef callMethods[] = {
    CALLDEF(dtw_basic, 10),
    CALLDEF(envelop, 2),
    CALLDEF(logGAK, 8),
    CALLDEF(pairs, 2),
    {NULL, NULL, 0}
};

void register_functions() {
    #define DTWCLUST_REGISTER(__FUN__) R_RegisterCCallable("dtwclust", #__FUN__, (DL_FUNC)__FUN__);
    DTWCLUST_REGISTER(dtw_basic)
    DTWCLUST_REGISTER(envelop)
    DTWCLUST_REGISTER(logGAK)
    DTWCLUST_REGISTER(pairs)
    #undef DTWCLUST_REGISTER
}

extern "C" void R_init_dtwclust(DllInfo* info) {
    register_functions();
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
