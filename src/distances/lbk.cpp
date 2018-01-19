#include "distances.h"

#include <cmath> // std::sqrt

#include <RcppArmadillo.h>

#include "../utils/utils.h" // kahan_sum

namespace dtwclust {

// thread-safe overload
double lbk_core(const double * const x, const int length, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const H)
{
    double lb = 0;
    for (int i = 0; i < length; i++) {
        if (x[i] > upper_envelope[i])
            H[i] = x[i] - upper_envelope[i];
        else if (x[i] < lower_envelope[i])
            H[i] = lower_envelope[i] - x[i];
        else
            H[i] = 0;
        if (p > 1) H[i] *= H[i];
    }
    lb = kahan_sum(H, length);
    if (p > 1) lb = std::sqrt(lb);
    return lb;
}

// non-thread-safe
double lbk_core(const Rcpp::NumericVector& x, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& H)
{
    return lbk_core(&x[0], x.length(), p, &lower_envelope[0], &upper_envelope[0], &H[0]);
}

// gateway
RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U)
{
    BEGIN_RCPP
    Rcpp::NumericVector x(X);
    Rcpp::NumericVector H(x.length());
    return Rcpp::wrap(lbk_core(x, Rcpp::as<int>(P), L, U, H));
    END_RCPP
}

} // namespace dtwclust
