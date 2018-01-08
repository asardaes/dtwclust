#include "distances++.h"

#include <cmath> // std::sqrt

#include <RcppArmadillo.h>

#include "../utils/utils++.h" // kahan_sum

namespace dtwclust {

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

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U) {
    BEGIN_RCPP
    Rcpp::NumericVector x(X);
    Rcpp::NumericVector H(x.length());
    return Rcpp::wrap(lbk_core(x, Rcpp::as<int>(P), L, U, H));
    END_RCPP
}

} // namespace dtwclust
