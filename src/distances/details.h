#ifndef DTWCLUST_DISTANCES_DETAILS_HPP_
#define DTWCLUST_DISTANCES_DETAILS_HPP_

#include <type_traits> // conditional

#include "../utils/SurrogateMatrix.h"
#include "../utils/utils.h" // id_t

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

// do not use volatile in x64, in x32 it's to avoid some comparison problems in dtw_basic which_min
typedef std::conditional<sizeof(void*) == 4, volatile double, double>::type dtwclust_tuple_t;

// dtw-basic.cpp
int backtrack_steps(const SurrogateMatrix<double>& lcm,
                    SurrogateMatrix<int>& index1,
                    SurrogateMatrix<int>& index2,
                    const std::size_t nx,
                    const std::size_t ny);
double dtw_basic_c(SurrogateMatrix<double>& lcm,
                   const SurrogateMatrix<const double>& x,
                   const SurrogateMatrix<const double>& y,
                   const int w,
                   const double norm,
                   const double step,
                   const bool backtrack);
double dtw_basic_par(SurrogateMatrix<double>& lcm,
                     const SurrogateMatrix<const double>& x,
                     const SurrogateMatrix<const double>& y,
                     const int window,
                     const double norm,
                     const double step,
                     const bool normalize);
double dtw_basic_par(SurrogateMatrix<double>& lcm,
                     const SurrogateMatrix<const double>& x,
                     const SurrogateMatrix<const double>& y,
                     const int window,
                     const double norm,
                     const double step,
                     const bool normalize,
                     SurrogateMatrix<int>& index1,
                     SurrogateMatrix<int>& index2,
                     int& path);

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
double logGAK_c(const SurrogateMatrix<const double>& seq1 ,
                const SurrogateMatrix<const double>& seq2,
                const double sigma,
                const id_t triangular,
                SurrogateMatrix<double>& logs);

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