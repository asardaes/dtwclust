#ifndef DTWCLUST_DISTMAT_LOOPS_HPP_
#define DTWCLUST_DISTMAT_LOOPS_HPP_

#define R_NO_REMAP
#include <Rinternals.h>
#undef R_NO_REMAP

namespace dtwclust {

extern "C" {
    SEXP distmat_loop(SEXP D, SEXP X, SEXP Y,
                      SEXP DIST, SEXP DIST_ARGS,
                      SEXP FILL_TYPE, SEXP MAT_TYPE, SEXP NUM_THREADS);

    SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS, SEXP NUM_THREADS);
}

} // namespace dtwclust

#endif // DTWCLUST_DISTMAT_LOOPS_HPP_
