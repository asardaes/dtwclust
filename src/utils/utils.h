#ifdef __cplusplus
extern "C" {
#endif

#ifndef DTWCLUST_UTILS_H_
#define DTWCLUST_UTILS_H_

#include <Rinternals.h>

SEXP pairs(SEXP L);

SEXP setnames_inplace(SEXP vec, SEXP names);

void Rflush();

// double to single index for matrices
int inline d2s(int const i, int const j, int const num_rows)
    __attribute__((always_inline));
int inline d2s(int const i, int const j, int const num_rows)
{ return i + j * num_rows; }

#endif // DTWCLUST_UTILS_H_

#ifdef __cplusplus
}
#endif
