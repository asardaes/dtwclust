#ifndef DTWCLUST_UTILS_HPP_
#define DTWCLUST_UTILS_HPP_

#include <cstddef> // std::size_t

namespace dtwclust {

typedef std::size_t id_t;

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

#define DTWCLUST_MIN_GRAIN 8

// envelope.cpp
void envelope_cpp(const double * const array, const int length, const unsigned int width,
                  double * const minvalues, double * const maxvalues);

// utils.cpp
void Rflush();
int get_grain(const int n, const int num_threads);
double kahan_sum(const double * const x, const int length);
void s2d(const int id, const int nrow, int& i, int& j);

} // namespace dtwclust

#endif // DTWCLUST_UTILS_HPP_
