#ifndef DTWCLUST_DISTANCES_HPP_
#define DTWCLUST_DISTANCES_HPP_

#define R_NO_REMAP
#include <Rinternals.h>
#undef R_NO_REMAP

namespace dtwclust {

extern "C" {
    SEXP dtw_basic(SEXP X, SEXP Y, SEXP WINDOW,
                   SEXP X_LEN, SEXP Y_LEN, SEXP NUM_VAR,
                   SEXP NORM, SEXP STEP, SEXP BACKTRACK, SEXP NORMALIZE, SEXP SQRT_DIST,
                   SEXP LCM);

    SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

    SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

    SEXP logGAK(SEXP X, SEXP Y, SEXP NX, SEXP NY, SEXP NUM_VAR,
                SEXP SIGMA, SEXP WINDOW, SEXP LOGS);

    SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV);
}

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_HPP_
