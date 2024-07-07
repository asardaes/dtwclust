#include "R-gateways.h"

#include <RcppArmadillo.h>

#include "../distances/calculators.h"
#include "distmat.h"
#include "fillers.h"

namespace dtwclust {

static auto distmat_factory = DistmatFactory();
static auto distance_calculator_factory = DistanceCalculatorFactory();
static auto distmat_filler_factory = DistmatFillerFactory();

extern "C" SEXP distmat_loop(SEXP D, SEXP X, SEXP Y,
                             SEXP DIST, SEXP DIST_ARGS,
                             SEXP FILL_TYPE, SEXP MAT_TYPE, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    auto distmat = distmat_factory.create(MAT_TYPE, D);
    auto dist_calculator = distance_calculator_factory.create(DIST, DIST_ARGS, X, Y);
    auto distmat_filler = distmat_filler_factory.create(FILL_TYPE, NUM_THREADS, distmat, dist_calculator);
    distmat_filler->fill();
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
