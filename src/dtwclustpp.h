#include <R.h>
#include <Rinternals.h>
#include <Rcpp.h>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

namespace dtwclust {

RcppExport SEXP envelope(SEXP series, SEXP window);

} // namespace dtwclust

#endif
