#include "distmat-fillers.h"

#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>

#include "concrete-fillers.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistmatFiller>
DistmatFillerFactory::create(const SEXP& FILL_TYPE,
                             const SEXP& NUM_THREADS,
                             std::shared_ptr<Distmat>& distmat,
                             const std::shared_ptr<DistanceCalculator>& dist_calculator)
{
    std::string type = Rcpp::as<std::string>(FILL_TYPE);
    if (type == "PAIRWISE") {
        return std::make_shared<PairwiseFiller>(distmat, dist_calculator, NUM_THREADS);
    }
    else if (type == "SYMMETRIC") {
        return std::make_shared<SymmetricFiller>(distmat, dist_calculator, NUM_THREADS);
    }
    else {
        return std::make_shared<PrimaryFiller>(distmat, dist_calculator, NUM_THREADS);
    }
}

} // namespace dtwclust
