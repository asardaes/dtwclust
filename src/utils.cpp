#include "dtwclust++.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* core logic to fill distmat in the different loops */
// =================================================================================================

void fill_distmat(DistmatFiller& distmat_filler,
                  const SEXP& D, const SEXP& X, const SEXP& Y,
                  const bool pairwise, const bool symmetric)
{
    distmat_filler.chooseFillStrategy(pairwise, symmetric);
    distmat_filler.fillDistmat(D, X, Y);
}

// =================================================================================================
/* for Rcpp::Rcout */
// =================================================================================================

void Rflush()
{
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
}

} // namespace dtwclust
