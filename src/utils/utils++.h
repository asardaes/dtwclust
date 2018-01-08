#ifndef DTWCLUST_UTILS_HPP_
#define DTWCLUST_UTILS_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* can be called from R */
// -------------------------------------------------------------------------------------------------

// envelope.cpp
RcppExport SEXP envelope(SEXP series, SEXP window);

// utils-cpp.cpp
RcppExport SEXP force_lb_symmetry(SEXP X);

// for sparse matrices in R
RcppExport SEXP SparseDistmatIndices__new(SEXP num_rows);
RcppExport SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions R */
// -------------------------------------------------------------------------------------------------

// utils-cpp.cpp
double kahan_sum(const Rcpp::NumericVector& x);

// envelope.cpp
void envelope_cpp(const Rcpp::NumericVector& array, const unsigned int width,
                  Rcpp::NumericVector& minvalues, Rcpp::NumericVector& maxvalues);

} // namespace dtwclust

#endif // DTWCLUST_UTILS_HPP_
