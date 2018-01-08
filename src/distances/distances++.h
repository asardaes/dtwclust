#ifndef DTWCLUST_DISTANCES_HPP_
#define DTWCLUST_DISTANCES_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP DISTMAT, SEXP MV);

double lbk_core(const Rcpp::NumericVector& x, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& H);

double lbi_core(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
                const unsigned int window_size, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& L2,
                Rcpp::NumericVector& U2,
                Rcpp::NumericVector& H,
                Rcpp::NumericVector& LB);

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_HPP_
