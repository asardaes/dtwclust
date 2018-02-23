#ifndef DTWCLUST_DISTANCES_HPP_
#define DTWCLUST_DISTANCES_HPP_

#include <RcppArmadillo.h>

#include "../utils/SurrogateMatrix.h"

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* can be called from R */
// -------------------------------------------------------------------------------------------------

RcppExport SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
                          SEXP m, SEXP n, SEXP num_var,
                          SEXP norm, SEXP step, SEXP backtrack, SEXP normalize,
                          SEXP distmat);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var,
                       SEXP sigma, SEXP window, SEXP logs);

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV);

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_HPP_
