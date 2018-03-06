#ifndef DTWCLUST_DISTANCES_DETAILS_HPP_
#define DTWCLUST_DISTANCES_DETAILS_HPP_

#include <type_traits> // conditional

#include "../utils/SurrogateMatrix.h"

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

// do not use volatile in x64, in x32 it's to avoid some comparison problems in dtw_basic which_min
typedef std::conditional<sizeof(void*) == 4, volatile double, double>::type dtwclust_tuple_t;

// dtw-basic.cpp
int backtrack_steps(double const * const D, int const nx, int const ny, int *index1, int *index2);
double dtw_basic_c(double * const D, dtwclust_tuple_t * const tuple,
                   double const * const x, double const * const y, int const w,
                   int const nx, int const ny, int const num_var,
                   double const norm, double const step,
                   int const backtrack);
double dtw_basic_par(double const * const x, double const * const y,
                     int const nx, int const ny, int const num_var,
                     int const window, double const norm, double const step, int const normalize,
                     double * const distmat);
double dtw_basic_par(double const * const x, double const * const y,
                     int const nx, int const ny, int const num_var,
                     int const window, double const norm, double const step, int const normalize,
                     double * const distmat,
                     int * const index1, int * const index2, int * const path);

// lbi.cpp
double lbi_core(const double * const x, const double * const y,
                const int length, const unsigned int window_size, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const L2, double * const U2, double * const H, double * const LB);

// lbk.cpp
double lbk_core(const double * const x, const int length, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const H);

// logGAK.cpp
double logGAK_c(double const *seq1 , double const *seq2,
                int const nX, int const nY, int const num_var,
                double const sigma, int const triangular,
                double *logs);
double logGAK_par(double const * const x, double const * const y,
                  int const nx, int const ny, int const num_var,
                  double const sigma, int const triangular,
                  double * const logs);

// soft-dtw.cpp
double soft_min(double a, double b, double c, const double gamma);
double squared_euclidean(const double * const x, const double * const y,
                         const int i, const int j,
                         const int x_nrows, const int y_nrows, const int ncols);
double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, SurrogateMatrix<double>& costmat);
double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, SurrogateMatrix<double>& costmat,
            SurrogateMatrix<double>& distmat);

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_DETAILS_HPP_
