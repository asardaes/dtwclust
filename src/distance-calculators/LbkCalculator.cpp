#include "distance-calculators.h"

#include <RcppArmadillo.h>

#include "../distances/distances++.h" // lbk_core

namespace dtwclust {

// =================================================================================================
/* lb_keogh distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbkCalculator::LbkCalculator(
    const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : DistanceCalculator(DIST_ARGS, X, Y)
{
    int length = Rcpp::as<int>((SEXP)dist_args_["len"]);
    p_ = Rcpp::as<int>((SEXP)dist_args_["p"]);
    lower_envelopes_ = dist_args_["lower.env"];
    upper_envelopes_ = dist_args_["upper.env"];
    H_ = Rcpp::NumericVector(length);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const Rcpp::NumericVector& x,
                                        const Rcpp::NumericVector& lower_envelope,
                                        const Rcpp::NumericVector& upper_envelope)
{
    return lbk_core(x, p_, lower_envelope, upper_envelope, H_);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const int i, const int j)
{
    // Y is ignored here, only the envelopes matter
    SEXP x = x_[i];
    SEXP lower_envelope = lower_envelopes_[j];
    SEXP upper_envelope = upper_envelopes_[j];
    return this->calculate(x, lower_envelope, upper_envelope);
}

} // namespace dtwclust
