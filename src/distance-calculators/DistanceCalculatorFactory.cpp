#include "distance-calculators.h"

#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>

#include "concrete-calculators.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::create(const SEXP& DIST, const SEXP& DIST_ARGS,
                                  const SEXP& X, const SEXP& Y)
{
    std::string dist = Rcpp::as<std::string>(DIST);
    return this->create(dist, DIST_ARGS, X, Y);
}

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::create(const std::string& dist, const SEXP& DIST_ARGS,
                                  const SEXP& X, const SEXP& Y)
{
    if (dist == "DTW_BASIC")
        return std::make_shared<DtwBasicCalculator>(DIST_ARGS, X, Y);
    else if (dist == "LBK")
        return std::make_shared<LbkCalculator>(DIST_ARGS, X, Y);
    else if (dist == "LBI")
        return std::make_shared<LbiCalculator>(DIST_ARGS, X, Y);
    else if (dist == "SDTW")
        return std::make_shared<SdtwCalculator>(DIST_ARGS, X, Y);
    else if (dist == "GAK")
        return std::make_shared<GakCalculator>(DIST_ARGS, X, Y);
    else if (dist == "SBD")
        return std::make_shared<SbdCalculator>(DIST_ARGS, X, Y);
    else
        Rcpp::stop("Unknown distance measure"); // nocov
}

} // namespace dtwclust
