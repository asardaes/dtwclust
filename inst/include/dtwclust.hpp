#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

namespace dtwclust {

SEXP dba(SEXP X, SEXP centroid,
         SEXP max_iter, SEXP delta, SEXP trace,
         SEXP multivariate, SEXP mv_ver, SEXP DOTS)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "dba");
    return fun(X, centroid, max_iter, delta, trace, multivariate, mv_ver, DOTS);
}

SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "dtw_lb");
    return fun(X, Y, D, MARGIN, DOTS);
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
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "soft_dtw");
    return fun(X, Y, GAMMA, COSTMAT, DISTMAT, MV);
}

SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "tadpole");
    return fun(X, K, DC, DTW_ARGS, LB, UB, TRACE, LIST);
}

} // namespace dtwclust

#endif
