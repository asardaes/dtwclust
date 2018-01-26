#ifndef DTWCLUST_UTILS_HPP_
#define DTWCLUST_UTILS_HPP_

#include <vector>

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

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

// utils.cpp
void Rflush();

// envelope.cpp
void envelope_cpp(const double * const array, const int length, const unsigned int width,
                  double * const minvalues, double * const maxvalues);

// utils.cpp
double kahan_sum(const double * const x, const int length);

// utils.cpp
void s2d(const int id, const int nrow, int& i, int& j);

// double to single index for matrices
int inline d2s(int const i, int const j, int const num_rows)
    __attribute__((always_inline));
int inline d2s(int const i, int const j, int const num_rows)
{ return i + j * num_rows; }

// grain parameter for multi-threading
int inline get_grain(const int n, const int num_threads)
    __attribute__((always_inline));
int inline get_grain(const int n, const int num_threads) {
    int grain = n / num_threads;
    // min grain defined here
    return (grain < 8) ? 8 : grain;
}

// for kahan sum (compensated sum)
class KahanSummer
{
public:
    KahanSummer(double * const x, const int nrows, const int ncols = 1);
    void reset();
    void add(const double value, const int i, const int j = 0);

private:
    double* const x_;
    int nrows_, ncols_;
    std::vector<double> c_, y_, t_;
};

} // namespace dtwclust

#endif // DTWCLUST_UTILS_HPP_
