#ifdef __cplusplus
extern "C" {
#endif

#ifndef DTWCLUST_DISTANCES_H_
#define DTWCLUST_DISTANCES_H_

#include <Rinternals.h>

SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP num_var,
               SEXP norm, SEXP step, SEXP backtrack,
               SEXP distmat);

double dtw_basic_par(double const * const x, double const * const y,
                     int const nx, int const ny, int const num_var,
                     int const window, double const norm, double const step,
                     double * const distmat);

SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var, SEXP sigma, SEXP window, SEXP logs);

#endif // DTWCLUST_DISTANCES_H_

#ifdef __cplusplus
}
#endif
