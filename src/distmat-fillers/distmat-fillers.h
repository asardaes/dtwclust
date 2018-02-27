#ifndef DTWCLUST_DISTMAT_FILLERS_HPP_
#define DTWCLUST_DISTMAT_FILLERS_HPP_

#include <memory> // *_ptr

#include <RcppArmadillo.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distmat/distmat.h"

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* abstract distmat filler */
// -------------------------------------------------------------------------------------------------
class DistmatFiller
{
public:
    virtual ~DistmatFiller() {};
    virtual void fill() const = 0;

protected:
    DistmatFiller(std::shared_ptr<Distmat>& distmat,
                  const std::shared_ptr<DistanceCalculator>& dist_calculator,
                  const SEXP& NUM_THREADS)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , num_threads_(Rcpp::as<int>(NUM_THREADS))
    { }

    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    std::shared_ptr<Distmat> distmat_;
    int num_threads_;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistmatFillerFactory
{
public:
    std::shared_ptr<DistmatFiller> create(
            const SEXP& FILL_TYPE,
            const SEXP& NUM_THREADS,
            std::shared_ptr<Distmat>& distmat,
            const std::shared_ptr<DistanceCalculator>& dist_calculator);
};

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

#endif // DTWCLUST_DISTMAT_FILLERS_HPP_
