#include <RcppArmadillo.h>

#include "distance-calculators.h"
#include "../distances/distances.h" // dtw_basic

namespace dtwclust {

// =================================================================================================
/* dtw_basic distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
DtwBasicDistanceCalculator::DtwBasicDistanceCalculator(
    const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : DistanceCalculator(DIST_ARGS, X, Y)
{
    window_ = dist_args_["window.size"];
    norm_ = dist_args_["norm"];
    step_ = dist_args_["step.pattern"];
    backtrack_ = dist_args_["backtrack"];
    gcm_ = dist_args_["gcm"];
    is_multivariate_ = Rcpp::as<bool>(dist_args_["is.multivariate"]);
    normalize_ = Rcpp::as<bool>(dist_args_["normalize"]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double DtwBasicDistanceCalculator::calculate(const SEXP& X, const SEXP& Y)
{
    bool backtrack = Rcpp::as<bool>(backtrack_);
    int x_len, y_len, num_vars;
    if (is_multivariate_) {
        Rcpp::NumericMatrix x(X), y(Y);
        x_len = x.nrow();
        y_len = y.nrow();
        num_vars = x.ncol();
    }
    else {
        Rcpp::NumericVector x(X), y(Y);
        x_len = x.length();
        y_len = y.length();
        num_vars = 1;
    }

    SEXP nx = PROTECT(Rcpp::wrap(x_len));
    SEXP ny = PROTECT(Rcpp::wrap(y_len));
    SEXP nv = PROTECT(Rcpp::wrap(num_vars));

    SEXP res = dtw_basic(X, Y, window_, nx, ny, nv, norm_, step_, backtrack_, gcm_);
    double distance;
    if (backtrack) { // nocov start
        Rcpp::List temp(res);
        distance = Rcpp::as<double>(temp["distance"]);
    } // nocov end
    else {
        distance = Rcpp::as<double>(res);
    }
    if (normalize_) {
        distance /= x_len + y_len;
    }
    UNPROTECT(3);
    return distance;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double DtwBasicDistanceCalculator::calculate(const int i, const int j)
{
    SEXP x = x_[i];
    SEXP y = y_[j];
    return this->calculate(x, y);
}

} // namespace dtwclust
