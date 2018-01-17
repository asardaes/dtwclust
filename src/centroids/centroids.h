#ifndef DTWCLUST_CENTROIDS_HPP_
#define DTWCLUST_CENTROIDS_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

RcppExport SEXP dba(SEXP X, SEXP CENT,
                    SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
                    SEXP MV, SEXP MV_VER, SEXP DOTS, SEXP NUM_THREADS);

RcppExport SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID,
                          SEXP GAMMA, SEXP WEIGHTS,
                          SEXP MV, SEXP NUM_THREADS);

} // namespace dtwclust

#endif // DTWCLUST_CENTROIDS_HPP_
