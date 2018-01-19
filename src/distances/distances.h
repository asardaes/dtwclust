#ifndef DTWCLUST_DISTANCES_HPP_
#define DTWCLUST_DISTANCES_HPP_

#include <RcppArmadillo.h>

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

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

double dtw_basic_par(double const * const x, double const * const y,
                     int const nx, int const ny, int const num_var,
                     int const window, double const norm, double const step, int const normalize,
                     double * const distmat, int const backtrack,
                     int * const index1, int * const index2, int * const path);

double lbi_core(const double * const x, const double * const y,
                const int length, const unsigned int window_size, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const L2, double * const U2, double * const H, double * const LB);

double lbk_core(const double * const x, const int length, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const H);

double logGAK_par(double const * const x, double const * const y,
                  int const nx, int const ny, int const num_var,
                  double const sigma, int const triangular,
                  double * const logs);

double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, double * const costmat,
            const bool save_norm, double * const distmat);

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_HPP_
