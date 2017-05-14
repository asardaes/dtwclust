#include <R.h>
#include <Rinternals.h>
#include <Rcpp.h>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

SEXP envelop(SEXP series, SEXP window)
{
    static SEXP(*fun)(SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("dtwclust", "envelop");
    return fun(series, window);
}

#endif
