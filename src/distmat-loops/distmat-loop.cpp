#include <RcppArmadillo.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distmat/distmat.h"
#include "../distmat-fillers/distmat-fillers.h"
#include "distmat-loops.h"

namespace dtwclust {

RcppExport SEXP distmat_loop(SEXP D, SEXP X, SEXP Y,
                             SEXP DIST, SEXP DIST_ARGS,
                             SEXP FILL_TYPE, SEXP MAT_TYPE, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    auto distmat = DistmatFactory().create(MAT_TYPE, D);
    auto dist_calculator = DistanceCalculatorFactory().create(DIST, DIST_ARGS, X, Y);
    auto distmat_filler = DistmatFillerFactory().create(FILL_TYPE, distmat, ENDPOINTS,
                                                        dist_calculator);
    distmat_filler->fill();
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
