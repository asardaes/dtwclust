#include "dtwclust++.h"

namespace dtwclust {

RcppExport SEXP lbi_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                         SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    DistmatFiller distmat_filler(BIGMAT, ENDPOINTS, Distance::LBI, DISTARGS);
    fill_distmat(distmat_filler, D, X, Y, Rcpp::as<bool>(PAIRWISE), Rcpp::as<bool>(SYMMETRIC));
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
