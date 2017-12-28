#include "dtwclust++.h"

namespace dtwclust {

RcppExport SEXP sbd_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                         SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    DistanceCalculatorFactory distcalc_factory;
    auto dist_calculator = distcalc_factory.createCalculator(Distance::SBD, DISTARGS);

    Distmat* distmat;
    if (Rcpp::as<bool>(BIGMAT)) {
        distmat = new BigmemoryDistmat(D);
    }
    else {
        distmat = new RDistmat(D);
    }

    DistmatFillerFactory distfill_factory;
    auto distmat_filler = distfill_factory.createFiller(
        Rcpp::as<bool>(PAIRWISE), Rcpp::as<bool>(SYMMETRIC),
        distmat, ENDPOINTS, dist_calculator
    );
    distmat_filler->fillDistmat(X, Y);

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
