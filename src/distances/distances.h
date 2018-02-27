#ifndef DTWCLUST_DISTANCES_HPP_
#define DTWCLUST_DISTANCES_HPP_

#define R_NO_REMAP
#include <Rinternals.h>

namespace dtwclust {

extern "C" {
    SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
                   SEXP m, SEXP n, SEXP num_var,
                   SEXP norm, SEXP step, SEXP backtrack, SEXP normalize,
                   SEXP distmat);

    SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

    SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

    SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var,
                SEXP sigma, SEXP window, SEXP logs);

    SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV);
}

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_HPP_
