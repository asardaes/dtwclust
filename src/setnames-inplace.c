#include <R.h>
#include <Rinternals.h>
#include "dtwclust.h"

SEXP setnames_inplace(SEXP vec, SEXP names) {
    setAttrib(vec, R_NamesSymbol, names);
    return R_NilValue;
}
