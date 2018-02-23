#ifndef DTWCLUST_R_UTILS_HPP_
#define DTWCLUST_R_UTILS_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* can be called from R */
// -------------------------------------------------------------------------------------------------

// for sparse matrices in R
RcppExport SEXP SparseDistmatIndices__new(SEXP num_rows);
RcppExport SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

// envelope.cpp
RcppExport SEXP envelope(SEXP series, SEXP window);

// utils.cpp
RcppExport SEXP force_lb_symmetry(SEXP X);

// utils.cpp
RcppExport SEXP pairs(SEXP L);

// utils.cpp
RcppExport SEXP setnames_inplace(SEXP vec, SEXP names);

// helper for semi-supervised dtwclust
RcppExport SEXP PairTracker__new(SEXP max_size);
RcppExport SEXP PairTracker__link(SEXP xptr, SEXP i, SEXP j, SEXP link_type);
RcppExport SEXP PairTracker__getUnseenPair(SEXP xptr);

} // namespace dtwclust

#endif // DTWCLUST_R_UTILS_HPP_
