#include <Rcpp.h>
#include <cmath>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* helper kahan_sum */
// =================================================================================================

double kahan_sum(const Rcpp::NumericVector& x) {
    double sum = 0, c = 0;
    for (double i : x) {
        double y = i - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

// =================================================================================================
/* LB_Keogh */
// =================================================================================================

SEXP lbk_cpp(const Rcpp::NumericVector& x, int p,
             const Rcpp::NumericVector& lower_envelope, const Rcpp::NumericVector& upper_envelope,
             Rcpp::NumericVector& H)
{
    double lb = 0;
    for (int i = 0; i < x.length(); i++) {
        if (x[i] > upper_envelope[i])
            H[i] = x[i] - upper_envelope[i];
        else if (x[i] < lower_envelope[i])
            H[i] = lower_envelope[i] - x[i];
        else
            H[i] = 0;

        if (p > 1) H[i] *= H[i];
    }

    lb = kahan_sum(H);
    if (p > 1) lb = std::sqrt(lb);
    return Rcpp::wrap(lb);
}

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U) {
BEGIN_RCPP
    Rcpp::NumericVector x(X), lower_envelope(L), upper_envelope(U);
    Rcpp::NumericVector H(x.length());
    return lbk_cpp(x, Rcpp::as<int>(P), lower_envelope, upper_envelope, H);
END_RCPP
}

// =================================================================================================
/* LB_Improved */
// =================================================================================================

SEXP lbi_cpp(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
             unsigned int window_size, int p,
             const Rcpp::NumericVector& lower_envelope, const Rcpp::NumericVector& upper_envelope,
             Rcpp::NumericVector& L2, Rcpp::NumericVector& U2, Rcpp::NumericVector&H)
{
    Rcpp::NumericVector LB(x.length());
    double lb = 0;
    for (int i = 0; i < x.length(); i++) {
        if (x[i] > upper_envelope[i]) {
            H[i] = upper_envelope[i];
            LB[i] = x[i] - upper_envelope[i];

        } else if (x[i] < lower_envelope[i]) {
            H[i] = lower_envelope[i];
            LB[i] = lower_envelope[i] - x[i];

        } else {
            H[i] = x[i];
            LB[i] = std::abs(x[i] - x[i]);
        }

        if (p > 1) LB[i] *= LB[i];
    }

    envelope_cpp(H, window_size * 2 + 1, L2, U2);

    double temp = 0;
    for (int i = 0; i < y.length(); i++) {
        if (y[i] > U2[i])
            temp = y[i] - U2[i];
        else if (y[i] < L2[i])
            temp = L2[i] - y[i];
        else
            temp = std::abs(y[i] - y[i]);

        if (p > 1) temp *= temp;
        LB[i] += temp;
    }

    lb = kahan_sum(LB);
    if (p > 1) lb = std::sqrt(lb);
    return Rcpp::wrap(lb);
}

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U) {
BEGIN_RCPP
    Rcpp::NumericVector x(X), y(Y), lower_envelope(L), upper_envelope(U);
    Rcpp::NumericVector L2(x.length()), U2(x.length()), H(x.length());
    return lbi_cpp(x, y, Rcpp::as<unsigned int>(WINDOW), Rcpp::as<int>(P),
                   lower_envelope, upper_envelope, L2, U2, H);
END_RCPP
}

} // namespace dtwclust
