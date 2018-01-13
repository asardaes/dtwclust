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
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

// utils-cpp.cpp
double kahan_sum(const double * const x, const int length);

// envelope.cpp
void envelope_cpp(const double * const array, const int length, const unsigned int width,
                  double * const minvalues, double * const maxvalues);

// utils-cpp.cpp
void s2d(const int id, const int nrow, int& i, int& j);

int inline get_grain(const int n, const int num_threads)
    __attribute__((always_inline));
int inline get_grain(const int n, const int num_threads) {
    int grain = n / num_threads;
    // min_grain defined here
    return (grain < 10) ? 10 : grain;
}

} // namespace dtwclust

#endif // DTWCLUST_UTILS_HPP_
