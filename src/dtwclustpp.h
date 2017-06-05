#include <Rcpp.h>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

namespace dtwclust {

RcppExport SEXP dba(SEXP X, SEXP centroid,
                    SEXP max_iter, SEXP delta, SEXP trace,
                    SEXP multivariate, SEXP mv_ver, SEXP DOTS);

RcppExport SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP DOTS);

RcppExport SEXP envelope(SEXP series, SEXP window);

void envelope_cpp(const Rcpp::NumericVector& array, unsigned int width,
                  Rcpp::NumericVector& minvalues, Rcpp::NumericVector& maxvalues);

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

} // namespace dtwclust

#endif
