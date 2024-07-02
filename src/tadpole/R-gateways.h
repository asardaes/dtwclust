#ifndef DTWCLUST_TADPOLE_HPP_
#define DTWCLUST_TADPOLE_HPP_

#include <algorithm> // stable_sort
#include <numeric> // iota
#include <utility> // swap
#include <vector>

#define R_NO_REMAP
#include <Rinternals.h>
#undef R_NO_REMAP

namespace dtwclust {

extern "C" SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST, SEXP NUM_THREADS);

} // namespace dtwclust

#endif // DTWCLUST_TADPOLE_HPP_
