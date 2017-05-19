#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>
#include <Rinternals.h>

#ifndef _DTWCLUST_H
#define _DTWCLUST_H

SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP dim,
               SEXP norm, SEXP step, SEXP backtrack,
               SEXP distmat);

SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP dim, SEXP sigma, SEXP window, SEXP logs);

SEXP pairs(SEXP L, SEXP lower);

SEXP setnames_inplace(SEXP vec, SEXP names);

#endif

#ifdef __cplusplus
}
#endif
