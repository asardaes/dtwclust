#ifndef DTWCLUST_UTILS_HPP_
#define DTWCLUST_UTILS_HPP_

#include <cstddef> // std::size_t

#include "SurrogateMatrix.h"

namespace dtwclust {

typedef std::size_t id_t;

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

#define DTWCLUST_MIN_GRAIN 8

// envelope.cpp
void envelope_cpp(const SurrogateMatrix<double>& array, const unsigned int width,
                  SurrogateMatrix<double>& minvalues, SurrogateMatrix<double>& maxvalues);

// utils.cpp
void Rflush();
int get_grain(const int n, const int num_threads);
double kahan_sum(const SurrogateMatrix<double>& x);
void s2d(const id_t id, const id_t nrow, id_t& i, id_t& j);

} // namespace dtwclust

#endif // DTWCLUST_UTILS_HPP_
