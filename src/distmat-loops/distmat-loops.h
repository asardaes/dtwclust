#ifndef DTWCLUST_DISTMAT_LOOPS_HPP_
#define DTWCLUST_DISTMAT_LOOPS_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

RcppExport SEXP distmat_loop(SEXP D, SEXP X, SEXP Y,
                             SEXP DIST, SEXP DIST_ARGS,
                             SEXP FILL_TYPE, SEXP MAT_TYPE, SEXP ENDPOINTS);

RcppExport SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS);

} // namespace dtwclust

#endif // DTWCLUST_DISTMAT_LOOPS_HPP_
