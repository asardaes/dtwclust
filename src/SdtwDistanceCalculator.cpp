#include "dtwclust++.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* soft-DTW distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SdtwDistanceCalculator::SdtwDistanceCalculator(const SEXP& DIST_ARGS)
    : DistanceCalculator(DIST_ARGS)
{
    gamma_ = dist_args_["gamma"];
    costmat_ = dist_args_["cm"];
    mv_ = dist_args_["mv"];
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double SdtwDistanceCalculator::calculateDistance(const SEXP& X, const SEXP& Y)
{
    // 'distmat' parameter is always NULL in here
    return Rcpp::as<double>(soft_dtw(X, Y, gamma_, costmat_, R_NilValue, mv_));
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SdtwDistanceCalculator::calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                                                     const int i, const int j)
{
    SEXP x = X[i];
    SEXP y = Y[j];
    return this->calculateDistance(x, y);
}

} // namespace dtwclust
