#include <R.h>
#include <Rinternals.h>
#include <Rcpp.h>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

namespace dtwclust {

SEXP envelope(SEXP series, SEXP window)
{
    static SEXP(*fun)(SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("dtwclust", "envelope");
    return fun(series, window);
}

} // namespace dtwclust

#endif
