#include "dtwclust++.h"
#include <memory> // make_shared

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::createCalculator(enum Distance distance, const SEXP& DIST_ARGS)
{
    switch (distance)
    {
    case Distance::DTW_BASIC:
        return std::make_shared<DtwBasicDistanceCalculator>(DIST_ARGS);
    case Distance::LBK:
        return std::make_shared<LbkDistanceCalculator>(DIST_ARGS);
    case Distance::LBI:
        return std::make_shared<LbiDistanceCalculator>(DIST_ARGS);
    case Distance::SDTW:
        return std::make_shared<SdtwDistanceCalculator>(DIST_ARGS);
    case Distance::GAK:
        return std::make_shared<GakDistanceCalculator>(DIST_ARGS);
    case Distance::SBD:
        return std::make_shared<SbdDistanceCalculator>(DIST_ARGS);
    default:
        Rcpp::stop("Unknown distance measure");
    }
}

} // namespace dtwclust
