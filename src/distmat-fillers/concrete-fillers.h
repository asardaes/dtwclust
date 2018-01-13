#ifndef DTWCLUST_CONCRETE_FILLERS_HPP_
#define DTWCLUST_CONCRETE_FILLERS_HPP_

#include <memory> // *_ptr

#include <RcppArmadillo.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distmat/distmat.h"
#include "distmat-fillers.h"

namespace dtwclust {

// =================================================================================================
/* DistmatFillers (concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* pairwise filler */
// -------------------------------------------------------------------------------------------------
class PairwiseFiller : public DistmatFiller
{
public:
    PairwiseFiller(std::shared_ptr<Distmat>& distmat,
                   const std::shared_ptr<DistanceCalculator>& dist_calculator,
                   const SEXP& NUM_THREADS);
    void fill() const override;
};

// -------------------------------------------------------------------------------------------------
/* symmetric filler */
// -------------------------------------------------------------------------------------------------
class SymmetricFiller : public DistmatFiller
{
public:
    SymmetricFiller(std::shared_ptr<Distmat>& distmat,
                    const std::shared_ptr<DistanceCalculator>& dist_calculator,
                    const SEXP& NUM_THREADS);
    void fill() const override;
};

// -------------------------------------------------------------------------------------------------
/* primary filler */
// -------------------------------------------------------------------------------------------------
class PrimaryFiller : public DistmatFiller
{
public:
    PrimaryFiller(std::shared_ptr<Distmat>& distmat,
                  const std::shared_ptr<DistanceCalculator>& dist_calculator,
                  const SEXP& NUM_THREADS);
    void fill() const override;
};

} // namespace dtwclust

#endif // DTWCLUST_CONCRETE_FILLERS_HPP_
