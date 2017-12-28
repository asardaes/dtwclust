#ifdef __cplusplus
extern "C" {
#endif

#ifndef _DTWCLUST_DISTANCES_H_
#define _DTWCLUST_DISTANCES_H_

#include <R.h>
#include <Rinternals.h>

SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP num_var,
               SEXP norm, SEXP step, SEXP backtrack,
               SEXP distmat);

SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var, SEXP sigma, SEXP window, SEXP logs);

#endif // _DTWCLUST_DISTANCES_H_

#ifdef __cplusplus
}
#endif
