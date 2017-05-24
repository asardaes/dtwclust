#include <R.h>
#include <Rinternals.h>
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

SEXP envelope(SEXP series, SEXP window)
{
    static SEXP(*fun)(SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("dtwclust", "envelope");
    return fun(series, window);
}

} // namespace dtwclust

#endif
