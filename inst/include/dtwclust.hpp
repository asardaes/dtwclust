#include <Rcpp.h>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

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

SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP DOTS)
{
    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("dtwclust", "dtw_lb");
    return fun(X, Y, D, DOTS);
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

} // namespace dtwclust

#endif
