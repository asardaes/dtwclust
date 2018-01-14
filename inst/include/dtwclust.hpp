#ifndef DTWCLUST_HPP_
#define DTWCLUST_HPP_

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

namespace dtwclust {

SEXP dba(SEXP X, SEXP CENT,
         SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
         SEXP MV, SEXP MV_VER, SEXP DOTS, SEXP NUM_THREADS)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "dba");
    return fun(X, CENT, MAX_ITER, DELTA, TRACE, MV, MV_VER, DOTS, NUM_THREADS);
}

SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS, SEXP NUM_THREADS)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "dtw_lb");
    return fun(X, Y, D, MARGIN, DOTS, NUM_THREADS);
}

SEXP envelope(SEXP series, SEXP window)
{
    static SEXP(*fun)(SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("dtwclust", "envelope");
    return fun(series, window);
}

SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "lbk");
    return fun(X, P, L, U);
}

SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "lbi");
    return fun(X, Y, WINDOW, P, L, U);
}

SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID, SEXP GAMMA, SEXP WEIGHTS, SEXP MV,
               SEXP COSTMAT, SEXP DISTMAT, SEXP EM)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "sdtw_cent");
    return fun(SERIES, CENTROID, GAMMA, WEIGHTS, MV, COSTMAT, DISTMAT, EM);
}

SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP DISTMAT, SEXP MV)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "soft_dtw");
    return fun(X, Y, GAMMA, COSTMAT, DISTMAT, MV);
}

SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST, SEXP NUM_THREADS)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "tadpole");
    return fun(X, K, DC, DTW_ARGS, LB, UB, TRACE, LIST, NUM_THREADS);
}

} // namespace dtwclust

#endif // DTWCLUST_HPP_
