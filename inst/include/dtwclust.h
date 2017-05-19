#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>
#include <Rinternals.h>

#ifndef _DTWCLUST_H
#define _DTWCLUST_H

SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP dim,
               SEXP norm, SEXP step, SEXP backtrack,
               SEXP distmat)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "dtw_basic");
    return fun(x, y, window, m, n, dim, norm, step, backtrack, distmat);
}

SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP dim, SEXP sigma, SEXP window, SEXP logs)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "logGAK");
    return fun(x, y, nx, ny, dim, sigma, window, logs);
}

SEXP pairs(SEXP L, SEXP lower)
{
    static SEXP(*fun)(SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("dtwclust", "pairs");
    return fun(L, lower);
}

SEXP setnames_inplace(SEXP vec, SEXP names)
{
    static SEXP(*fun)(SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("dtwclust", "setnames_inplace");
    return fun(vec, names);
}

#endif

#ifdef __cplusplus
}
#endif
