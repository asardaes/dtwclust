#include "concrete-distmat.h"

#include <RcppArmadillo.h>
#include <RcppParallel.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
RDistmat::RDistmat(const SEXP& D)
    : distmat_(RcppParallel::RMatrix<double>(Rcpp::NumericMatrix(D)))
{ }

// -------------------------------------------------------------------------------------------------
/* operator() for assignment */
// -------------------------------------------------------------------------------------------------
double& RDistmat::operator() (const int i, const int j)
{ return distmat_(i,j); }

} // namespace dtwclust
