#ifndef DTWCLUST_R_UTILS_HPP_
#define DTWCLUST_R_UTILS_HPP_

#define R_NO_REMAP
#include <Rinternals.h>

namespace dtwclust {

extern "C" {
    SEXP SparseDistmatIndices__new(SEXP num_rows);
    SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

    SEXP envelope(SEXP series, SEXP window);
    SEXP force_lb_symmetry(SEXP X);
    SEXP pairs(SEXP L);
    SEXP setnames_inplace(SEXP vec, SEXP names);

    SEXP PairTracker__new(SEXP max_size);
    SEXP PairTracker__link(SEXP xptr, SEXP i, SEXP j, SEXP link_type);
    SEXP PairTracker__getUnseenPair(SEXP xptr);
}

} // namespace dtwclust

#endif // DTWCLUST_R_UTILS_HPP_
