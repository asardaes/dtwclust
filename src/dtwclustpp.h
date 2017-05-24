#include <R.h>
#include <Rinternals.h>
#include <Rcpp.h>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

namespace dtwclust {

RcppExport SEXP dba(SEXP X, SEXP centroid,
                    SEXP max_iter, SEXP delta, SEXP trace,
                    SEXP multivariate, SEXP mv_ver, SEXP DOTS);

RcppExport SEXP envelope(SEXP series, SEXP window);

} // namespace dtwclust

#endif
