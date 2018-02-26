#ifndef DTWCLUST_R_UTILS_HPP_
#define DTWCLUST_R_UTILS_HPP_

#define R_NO_REMAP
#include <Rinternals.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* can be called from R */
// -------------------------------------------------------------------------------------------------

extern "C" {
    // for sparse matrices in R
    SEXP SparseDistmatIndices__new(SEXP num_rows);
    SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

    // envelope.cpp
    SEXP envelope(SEXP series, SEXP window);

    // utils.cpp
    SEXP force_lb_symmetry(SEXP X);

    // utils.cpp
    SEXP pairs(SEXP L);

    // utils.cpp
    SEXP setnames_inplace(SEXP vec, SEXP names);

    // helper for semi-supervised dtwclust
    SEXP PairTracker__new(SEXP max_size);
    SEXP PairTracker__link(SEXP xptr, SEXP i, SEXP j, SEXP link_type);
    SEXP PairTracker__getUnseenPair(SEXP xptr);
}

} // namespace dtwclust

#endif // DTWCLUST_R_UTILS_HPP_
