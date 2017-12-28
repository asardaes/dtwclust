#include "dtwclust++.h"
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

namespace dtwclust {

// =================================================================================================
/* R matrix */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
RDistmat::RDistmat(const SEXP& D)
    : distmat_(Rcpp::NumericMatrix(D))
{ }

// -------------------------------------------------------------------------------------------------
/* operator() for assignment */
// -------------------------------------------------------------------------------------------------
double& RDistmat::operator() (const int i, const int j)
{
    return distmat_(i,j);
}

// =================================================================================================
/* big matrix */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor
 *   distmat_ **must** be initialized in the initializer list, becuase MatrixAccessor<double> has no
 *   default constructor (I think)
 */
// -------------------------------------------------------------------------------------------------
BigmemoryDistmat::BigmemoryDistmat(const SEXP& D)
    : distmat_(MatrixAccessor<double>(*(Rcpp::XPtr<BigMatrix>(D))))
{ }

// -------------------------------------------------------------------------------------------------
/* operator() for assignment */
// -------------------------------------------------------------------------------------------------
double& BigmemoryDistmat::operator() (const int i, const int j)
{
    // bigmemory operator[][] is backwards
    return distmat_[j][i];
}

} // namespace dtwclust
