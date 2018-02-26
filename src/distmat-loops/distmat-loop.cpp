#include "distmat-loops.h"

#include <RcppArmadillo.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distmat/distmat.h"
#include "../distmat-fillers/distmat-fillers.h"

namespace dtwclust {

extern "C" SEXP distmat_loop(SEXP D, SEXP X, SEXP Y,
                             SEXP DIST, SEXP DIST_ARGS,
                             SEXP FILL_TYPE, SEXP MAT_TYPE, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    auto distmat = DistmatFactory().create(MAT_TYPE, D);
    auto dist_calculator = DistanceCalculatorFactory()
        .create(DIST, DIST_ARGS, X, Y);
    auto distmat_filler = DistmatFillerFactory()
        .create(FILL_TYPE, NUM_THREADS, distmat, dist_calculator);
    distmat_filler->fill();
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
