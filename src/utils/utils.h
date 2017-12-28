#ifdef __cplusplus
extern "C" {
#endif

#ifndef _DTWCLUST_UTILS_H_
#define _DTWCLUST_UTILS_H_

#include <R.h>
#include <Rinternals.h>

SEXP pairs(SEXP L);

SEXP setnames_inplace(SEXP vec, SEXP names);

void Rflush();

#endif // _DTWCLUST_UTILS_H_

#ifdef __cplusplus
}
#endif
