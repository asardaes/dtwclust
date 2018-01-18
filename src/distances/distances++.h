#ifndef DTWCLUST_DISTANCES_HPP_
#define DTWCLUST_DISTANCES_HPP_

#include <RcppArmadillo.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* can be called from R */
// -------------------------------------------------------------------------------------------------

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV);

// -------------------------------------------------------------------------------------------------
/* called by other C++ functions */
// -------------------------------------------------------------------------------------------------

double lbk_core(const double * const x, const int length, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const H);

double lbi_core(const double * const x, const double * const y,
                const int length, const unsigned int window_size, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const L2, double * const U2, double * const H, double * const LB);

double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, double * const costmat,
            const bool save_norm, double * const distmat);

} // namespace dtwclust

#endif // DTWCLUST_DISTANCES_HPP_
