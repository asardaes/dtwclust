#ifdef __cplusplus
extern "C" {
#endif

#ifndef DTWCLUST_DISTANCES_H_
#define DTWCLUST_DISTANCES_H_

#include <Rinternals.h>

// -------------------------------------------------------------------------------------------------
/* can be called from R */
// -------------------------------------------------------------------------------------------------

SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP num_var,
               SEXP norm, SEXP step, SEXP backtrack, SEXP normalize,
               SEXP distmat);

SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var, SEXP sigma, SEXP window, SEXP logs);

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

double dtw_basic_par(double const * const x, double const * const y,
                     int const nx, int const ny, int const num_var,
                     int const window, double const norm, double const step, int const normalize,
                     double * const distmat);

double logGAK_par(double const * const x, double const * const y,
                  int const nx, int const ny, int const num_var,
                  double const sigma, int const triangular,
                  double * const logs);

#endif // DTWCLUST_DISTANCES_H_

#ifdef __cplusplus
}
#endif
