#include <Rcpp.h>
#include <cmath>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* helper kahan_sum */
// =================================================================================================

double kahan_sum(const Rcpp::NumericVector& x) {
    double sum = 0, c = 0;
    for (const double& i : x) {
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

double lbk_core(const Rcpp::NumericVector& x, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
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
    return lb;
}

double lbk_cpp(const Rcpp::NumericVector& x, const int p,
               const Rcpp::NumericVector& lower_envelope, const Rcpp::NumericVector& upper_envelope)
{
    Rcpp::NumericVector H(x.length());
    return lbk_core(x, p, lower_envelope, upper_envelope, H);
}

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U) {
    BEGIN_RCPP
    return Rcpp::wrap(lbk_cpp(X, Rcpp::as<int>(P), L, U));
    END_RCPP
}

// =================================================================================================
/* LB_Improved */
// =================================================================================================

double lbi_core(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
                const unsigned int window_size, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& L2,
                Rcpp::NumericVector& U2,
                Rcpp::NumericVector& H,
                Rcpp::NumericVector& LB)
{
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
            LB[i] = 0;
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
            temp = 0;

        if (p > 1) temp *= temp;
        LB[i] += temp;
    }

    lb = kahan_sum(LB);
    if (p > 1) lb = std::sqrt(lb);
    return lb;
}

double lbi_cpp(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
               const unsigned int window_size, const int p,
               const Rcpp::NumericVector& lower_envelope, const Rcpp::NumericVector& upper_envelope)
{
    Rcpp::NumericVector L2(x.length()), U2(x.length()), H(x.length());
    Rcpp::NumericVector LB(x.length());
    return lbi_core(x, y, window_size, p, lower_envelope, upper_envelope, L2, U2, H, LB);
}

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U) {
    BEGIN_RCPP
    return Rcpp::wrap(lbi_cpp(X, Y, Rcpp::as<unsigned int>(WINDOW), Rcpp::as<int>(P), L, U));
    END_RCPP
}

// =================================================================================================
/* Force symmetry helper */
// =================================================================================================

RcppExport SEXP force_lb_symmetry(SEXP X) {
    BEGIN_RCPP
    Rcpp::NumericMatrix matrix(X);
    for (int i = 1; i < matrix.nrow(); i++) {
        R_CheckUserInterrupt();
        for (int j = 0; j < i; j++) {
            double lb1 = matrix(i,j);
            double lb2 = matrix(j,i);
            if (lb1 > lb2)
                matrix(j,i) = lb1;
            else
                matrix(i,j) = lb2;
        }
    }
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
