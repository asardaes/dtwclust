#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* lb_improved distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbiDistanceCalculator::LbiDistanceCalculator(const SEXP& DIST_ARGS)
    : DistanceCalculator(DIST_ARGS)
{
    int length = Rcpp::as<int>((SEXP)dist_args_["len"]);
    p_ = Rcpp::as<int>((SEXP)dist_args_["p"]);
    window_size_ = Rcpp::as<unsigned int>((SEXP)dist_args_["window.size"]);
    lower_envelopes_ = dist_args_["lower.env"];
    upper_envelopes_ = dist_args_["upper.env"];
    H_ = Rcpp::NumericVector(length);
    L2_ = Rcpp::NumericVector(length);
    U2_ = Rcpp::NumericVector(length);
    LB_ = Rcpp::NumericVector(length);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbiDistanceCalculator::calculateDistance(const Rcpp::NumericVector& x,
                                                const Rcpp::NumericVector& y,
                                                const Rcpp::NumericVector& lower_envelope,
                                                const Rcpp::NumericVector& upper_envelope)
{
    return lbi_core(x, y, window_size_, p_, lower_envelope, upper_envelope,
                    L2_, U2_, H_, LB_);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double LbiDistanceCalculator::calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                                                const int i, const int j)
{
    SEXP x = X[i];
    SEXP y = Y[j];
    SEXP lower_envelope = lower_envelopes_[j];
    SEXP upper_envelope = upper_envelopes_[j];
    return this->calculateDistance(x, y, lower_envelope, upper_envelope);
}

} // namespace dtwclust
