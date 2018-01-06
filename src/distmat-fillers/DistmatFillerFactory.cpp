#include "distmat-fillers.h"

#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistmatFiller>
DistmatFillerFactory::create(const SEXP& FILL_TYPE,
                             std::shared_ptr<Distmat>& distmat, const SEXP& ENDPOINTS,
                             const std::shared_ptr<DistanceCalculator>& dist_calculator)
{
    string type = Rcpp::as<string>(FILL_TYPE);
    if (type == "PAIRWISE") {
        return std::make_shared<PairwiseDistmatFiller>(distmat, ENDPOINTS, dist_calculator);
    }
    else if (type == "SYMMETRIC") {
        return std::make_shared<SymmetricDistmatFiller>(distmat, ENDPOINTS, dist_calculator);
    }
    else {
        return std::make_shared<GeneralDistmatFiller>(distmat, ENDPOINTS, dist_calculator);
    }
}

} // namespace dtwclust
