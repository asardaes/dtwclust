#include "distances.h"

#include <cmath> // std::sqrt

#include <RcppArmadillo.h>

#include "../utils/utils.h" // envelope_cpp, kahan_sum

namespace dtwclust {

// thread-safe overload
double lbi_core(const double * const x, const double * const y,
                const int length, const unsigned int window_size, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const L2, double * const U2, double * const H, double * const LB)
{
    double lb = 0;
    for (int i = 0; i < length; i++) {
        if (x[i] > upper_envelope[i]) {
            H[i] = upper_envelope[i];
            LB[i] = x[i] - upper_envelope[i];
        }
        else if (x[i] < lower_envelope[i]) {
            H[i] = lower_envelope[i];
            LB[i] = lower_envelope[i] - x[i];
        }
        else {
            H[i] = x[i];
            LB[i] = 0;
        }
        if (p > 1) LB[i] *= LB[i];
    }
    envelope_cpp(H, length, window_size * 2 + 1, L2, U2);
    double temp = 0;
    for (int i = 0; i < length; i++) {
        if (y[i] > U2[i])
            temp = y[i] - U2[i];
        else if (y[i] < L2[i])
            temp = L2[i] - y[i];
        else
            temp = 0;
        if (p > 1) temp *= temp;
        LB[i] += temp;
    }
    lb = kahan_sum(LB, length);
    if (p > 1) lb = std::sqrt(lb);
    return lb;
}

// non-thread-safe
double lbi_core(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
                const unsigned int window_size, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& L2,
                Rcpp::NumericVector& U2,
                Rcpp::NumericVector& H,
                Rcpp::NumericVector& LB)
{
    return lbi_core(&x[0], &y[0], x.length(), window_size, p,
                    &lower_envelope[0], &upper_envelope[0],
                    &L2[0], &U2[0], &H[0], &LB[0]);
}

// gateway
RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U)
{
    BEGIN_RCPP
    Rcpp::NumericVector x(X);
    Rcpp::NumericVector L2(x.length()), U2(x.length()), H(x.length());
    Rcpp::NumericVector LB(x.length());
    return Rcpp::wrap(lbi_core(x, Y, Rcpp::as<unsigned int>(WINDOW), Rcpp::as<int>(P),
                               L, U, L2, U2, H, LB));
    END_RCPP
}

} // namespace dtwclust
