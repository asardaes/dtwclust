#include "dtwclust++.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* dtw_basic distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
DtwBasicDistanceCalculator::DtwBasicDistanceCalculator(const SEXP& DIST_ARGS)
    : DistanceCalculator(DIST_ARGS)
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
double DtwBasicDistanceCalculator::calculateDistance(const SEXP& X, const SEXP& Y)
{
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

    double distance = Rcpp::as<double>(
        dtw_basic(X, Y, window_, nx, ny, nv, norm_, step_, backtrack_, gcm_)
    );
    if (normalize_) {
        distance /= x_len + y_len;
    }
    UNPROTECT(3);
    return distance;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double DtwBasicDistanceCalculator::calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                                                     const int i, const int j)
{
    SEXP x = X[i];
    SEXP y = Y[j];
    return this->calculateDistance(x, y);
}

} // namespace dtwclust
