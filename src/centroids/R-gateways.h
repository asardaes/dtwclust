#ifndef DTWCLUST_CENTROIDS_HPP_
#define DTWCLUST_CENTROIDS_HPP_

#define R_NO_REMAP
#include <Rinternals.h>

namespace dtwclust {

extern "C" {
    SEXP dba(SEXP X, SEXP CENT,
             SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
             SEXP MV, SEXP MV_VER, SEXP DOTS, SEXP NUM_THREADS);

    SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID,
                             SEXP GAMMA, SEXP WEIGHTS,
                             SEXP MV, SEXP NUM_THREADS);
}

} // namespace dtwclust

#endif // DTWCLUST_CENTROIDS_HPP_
