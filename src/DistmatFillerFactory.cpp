#include "dtwclust++.h"
#include <memory> // make_shared

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistmatFiller>
DistmatFillerFactory::createFiller(const bool pairwise, const bool symmetric,
                                   Distmat* distmat, const SEXP& ENDPOINTS,
                                   const std::shared_ptr<DistanceCalculator>& dist_calculator)
{
    if (pairwise) {
        return std::make_shared<PairwiseDistmatFiller>(distmat, ENDPOINTS, dist_calculator);
    }
    else if (symmetric) {
        return std::make_shared<SymmetricDistmatFiller>(distmat, ENDPOINTS, dist_calculator);
    }
    else {
        return std::make_shared<GeneralDistmatFiller>(distmat, ENDPOINTS, dist_calculator);
    }
}

} // namespace dtwclust
