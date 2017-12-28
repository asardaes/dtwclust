#ifndef DTWCLUST_CENTROIDS_HPP_
#define DTWCLUST_CENTROIDS_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

RcppExport SEXP dba(SEXP X, SEXP centroid,
                    SEXP max_iter, SEXP delta, SEXP trace,
                    SEXP multivariate, SEXP mv_ver, SEXP DOTS);

RcppExport SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID, SEXP GAMMA, SEXP WEIGHTS, SEXP MV,
                          SEXP COSTMAT, SEXP DISTMAT, SEXP EM);

} // namespace dtwclust

#endif // DTWCLUST_CENTROIDS_HPP_
