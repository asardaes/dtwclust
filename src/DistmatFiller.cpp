#include "dtwclust++.h"
#include <memory> // *_ptr

namespace dtwclust {

// constructor
DistmatFiller::DistmatFiller(const SEXP& IS_BIGMAT, const SEXP& ENDPOINTS,
                             enum Distance distance, const SEXP& DIST_ARGS)
    : is_bigmat_(Rcpp::as<bool>(IS_BIGMAT))
    , ENDPOINTS_(ENDPOINTS)
{
    DistanceCalculatorFactory factory;
    dist_calculator_ = factory.createCalculator(distance, DIST_ARGS);
}

// chooseFillStrategy
void DistmatFiller::chooseFillStrategy(const bool pairwise, const bool symmetric)
{
    if (pairwise) {
        fill_strategy_ = std::unique_ptr<PairwiseDistmatFill>(new PairwiseDistmatFill());
    }
    else if (symmetric) {
        fill_strategy_ = std::unique_ptr<SymmetricDistmatFill>(new SymmetricDistmatFill());
    }
    else {
        fill_strategy_ = std::unique_ptr<GeneralDistmatFill>(new GeneralDistmatFill());
    }
}

// fillDistmat
void DistmatFiller::fillDistmat(const SEXP& D, const SEXP& X, const SEXP& Y) const
{
    if (fill_strategy_ != nullptr) {
        fill_strategy_->fillDistmat(D, X, Y, dist_calculator_, ENDPOINTS_, is_bigmat_);
    }
    else { // nocov start
        Rcpp::stop("Attempting to fill a distance matrix without choosing strategy");
    } // nocov end
}

} // namespace dtwclust
