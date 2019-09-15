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
double dtw_basic(SurrogateMatrix<double>& lcm,
                 const SurrogateMatrix<const double>& x,
                 const SurrogateMatrix<const double>& y,
                 const int window,
                 const double norm,
                 const double step,
                 const bool normalize,
                 const bool sqrt_dist);

double dtw_basic(SurrogateMatrix<double>& lcm,
                 const SurrogateMatrix<const double>& x,
                 const SurrogateMatrix<const double>& y,
                 const int window,
                 const double norm,
                 const double step,
                 const bool normalize,
                 const bool sqrt_dist,
                 SurrogateMatrix<int>& index1,
                 SurrogateMatrix<int>& index2,
                 int& path);

// lbi.cpp
double lbi_core(const SurrogateMatrix<const double>& x,
                const SurrogateMatrix<const double>& y,
                const unsigned int window_size,
                const int p,
                const SurrogateMatrix<const double>& lower_envelope,
                const SurrogateMatrix<const double>& upper_envelope,
                SurrogateMatrix<double>& L2,
                SurrogateMatrix<double>& U2,
                SurrogateMatrix<double>& H,
                SurrogateMatrix<double>& LB);

// lbk.cpp
double lbk_core(const SurrogateMatrix<const double>& x,
                const int p,
                const SurrogateMatrix<const double>& lower_envelope,
                const SurrogateMatrix<const double>& upper_envelope,
                SurrogateMatrix<double>& H);

// logGAK.cpp
double logGAK_c(const SurrogateMatrix<const double>& seq1 ,
                const SurrogateMatrix<const double>& seq2,
                const double sigma,
                const id_t triangular,
                SurrogateMatrix<double>& logs);

// soft-dtw.cpp
double sdtw(const SurrogateMatrix<const double>& x, const SurrogateMatrix<const double>& y,
            const double gamma, SurrogateMatrix<double>& costmat);
double sdtw(const SurrogateMatrix<const double>& x, const SurrogateMatrix<const double>& y,
            const double gamma, SurrogateMatrix<double>& costmat, SurrogateMatrix<double>& distmat);

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_DETAILS_HPP_
