#include "dtwclust++.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* dtw_basic distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
GakDistanceCalculator::GakDistanceCalculator(const SEXP& DIST_ARGS)
    : DistanceCalculator(DIST_ARGS)
{
    sigma_ = dist_args_["sigma"];
    window_ = dist_args_["window.size"];
    logs_ = dist_args_["logs"];
    is_multivariate_ = Rcpp::as<bool>(dist_args_["is.multivariate"]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double GakDistanceCalculator::calculate(const SEXP& X, const SEXP& Y)
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
        logGAK(X, Y, nx, ny, nv, sigma_, window_, logs_)
    );
    UNPROTECT(3);
    return distance;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double GakDistanceCalculator::calculate(const Rcpp::List& X, const Rcpp::List& Y,
                                        const int i, const int j)
{
    SEXP x = X[i];
    SEXP y = Y[j];
    return this->calculate(x, y);
}

} // namespace dtwclust
