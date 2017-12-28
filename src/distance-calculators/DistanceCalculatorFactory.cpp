#include <RcppArmadillo.h>
#include <memory> // make_shared
#include <string>

#include "distance-calculators.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::create(const SEXP& DIST, const SEXP& DIST_ARGS)
{
    std::string dist = Rcpp::as<std::string>(DIST);
    return this->create(dist, DIST_ARGS);
}

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::create(const std::string& dist, const SEXP& DIST_ARGS)
{
    if (dist == "DTW_BASIC")
        return std::make_shared<DtwBasicDistanceCalculator>(DIST_ARGS);
    else if (dist == "LBK")
        return std::make_shared<LbkDistanceCalculator>(DIST_ARGS);
    else if (dist == "LBI")
        return std::make_shared<LbiDistanceCalculator>(DIST_ARGS);
    else if (dist == "SDTW")
        return std::make_shared<SdtwDistanceCalculator>(DIST_ARGS);
    else if (dist == "GAK")
        return std::make_shared<GakDistanceCalculator>(DIST_ARGS);
    else if (dist == "SBD")
        return std::make_shared<SbdDistanceCalculator>(DIST_ARGS);
    else
        Rcpp::stop("Unknown distance measure"); // nocov
}

} // namespace dtwclust
