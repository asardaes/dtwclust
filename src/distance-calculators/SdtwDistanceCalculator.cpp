#include "distance-calculators.h"

#include <RcppArmadillo.h>

#include "../distances/distances++.h" // soft_dtw

namespace dtwclust {

// =================================================================================================
/* soft-DTW distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SdtwDistanceCalculator::SdtwDistanceCalculator(
    const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : DistanceCalculator(DIST_ARGS, X, Y)
{
    gamma_ = dist_args_["gamma"];
    costmat_ = dist_args_["cm"];
    mv_ = dist_args_["mv"];
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double SdtwDistanceCalculator::calculate(const SEXP& X, const SEXP& Y)
{
    // 'distmat' parameter is always NULL in here
    return Rcpp::as<double>(soft_dtw(X, Y, gamma_, costmat_, R_NilValue, mv_));
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SdtwDistanceCalculator::calculate(const int i, const int j)
{
    SEXP x = x_[i];
    SEXP y = y_[j];
    return this->calculate(x, y);
}

} // namespace dtwclust
