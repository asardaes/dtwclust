#ifndef DTWCLUST_DISTMAT_HPP_
#define DTWCLUST_DISTMAT_HPP_

#include <memory> // *_ptr

#include <RcppArmadillo.h>
#include <RcppParallel.h>

namespace dtwclust {

// =================================================================================================
/* Distmat (base + factory) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distmat */
// -------------------------------------------------------------------------------------------------
class Distmat
{
public:
    virtual ~Distmat() {};
    virtual double& operator() (const int i, const int j) = 0;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistmatFactory
{
public:
    std::shared_ptr<Distmat> create(const SEXP& MAT_TYPE, const SEXP& D);
};

} // namespace dtwclust

#endif // DTWCLUST_DISTMAT_HPP_
