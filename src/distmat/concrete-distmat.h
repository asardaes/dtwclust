#ifndef DTWCLUST_CONCRETE_DISTMAT_HPP_
#define DTWCLUST_CONCRETE_DISTMAT_HPP_

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "distmat.h"

namespace dtwclust {

// =================================================================================================
/* Distmat (concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* R matrix distmat (thread-safe) */
// -------------------------------------------------------------------------------------------------
class RDistmat : public Distmat
{
public:
    RDistmat(const SEXP& D);
    double& operator() (const int i, const int j) override;

private:
    RcppParallel::RMatrix<double> distmat_;
};

} // namespace dtwclust

#endif // DTWCLUST_CONCRETE_DISTMAT_HPP_
