#ifndef DTWCLUST_MISC_HPP_
#define DTWCLUST_MISC_HPP_

#include <R.h>
#include <Rinternals.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* for sparse matrices in R */
// -------------------------------------------------------------------------------------------------

RcppExport SEXP SparseDistmatIndices__new(SEXP num_rows);
RcppExport SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

} // namespace dtwclust

#endif // DTWCLUST_MISC_HPP_
