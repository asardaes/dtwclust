#ifndef _DTWCLUST_DISTMAT_FILLERS_HPP_
#define _DTWCLUST_DISTMAT_FILLERS_HPP_

#include <RcppArmadillo.h>
#include <memory> // *_ptr

#include "../distance-calculators/distance-calculators.h"
#include "../distmat/distmat.h"

namespace dtwclust {

// =================================================================================================
/* DistmatFillers (base + factory + concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distmat filler */
// -------------------------------------------------------------------------------------------------
class DistmatFiller
{
public:
    virtual ~DistmatFiller() {};
    virtual void fill(const Rcpp::List& X, const Rcpp::List& Y) const = 0;

protected:
    DistmatFiller(std::shared_ptr<Distmat>& distmat, const SEXP& ENDPOINTS,
                  const std::shared_ptr<DistanceCalculator>& dist_calculator)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , endpoints_(ENDPOINTS)
    { }

    std::shared_ptr<DistanceCalculator> dist_calculator_;
    std::shared_ptr<Distmat> distmat_;
    SEXP endpoints_;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistmatFillerFactory
{
public:
    std::shared_ptr<DistmatFiller>
    create(const SEXP& FILL_TYPE, std::shared_ptr<Distmat>& distmat, const SEXP& ENDPOINTS,
           const std::shared_ptr<DistanceCalculator>& dist_calculator);
};

// -------------------------------------------------------------------------------------------------
/* pairwise distmat filler */
// -------------------------------------------------------------------------------------------------
class PairwiseDistmatFiller : public DistmatFiller
{
public:
    PairwiseDistmatFiller(std::shared_ptr<Distmat>& distmat, const SEXP& ENDPOINTS,
                          const std::shared_ptr<DistanceCalculator>& dist_calculator)
        : DistmatFiller(distmat, ENDPOINTS, dist_calculator)
    { }

    void fill(const Rcpp::List& X, const Rcpp::List& Y) const override;
};

// -------------------------------------------------------------------------------------------------
/* symmetric distmat filler */
// -------------------------------------------------------------------------------------------------
class SymmetricDistmatFiller : public DistmatFiller
{
public:
    SymmetricDistmatFiller(std::shared_ptr<Distmat>& distmat, const SEXP& ENDPOINTS,
                           const std::shared_ptr<DistanceCalculator>& dist_calculator)
        : DistmatFiller(distmat, ENDPOINTS, dist_calculator)
    { }

    void fill(const Rcpp::List& X, const Rcpp::List& Y) const override;
};

// -------------------------------------------------------------------------------------------------
/* general distmat filler */
// -------------------------------------------------------------------------------------------------
class GeneralDistmatFiller : public DistmatFiller
{
public:
    GeneralDistmatFiller(std::shared_ptr<Distmat>& distmat, const SEXP& ENDPOINTS,
                         const std::shared_ptr<DistanceCalculator>& dist_calculator)
        : DistmatFiller(distmat, ENDPOINTS, dist_calculator)
    { }

    void fill(const Rcpp::List& X, const Rcpp::List& Y) const override;
};

} // namespace dtwclust

#endif // _DTWCLUST_DISTMAT_FILLERS_HPP_
