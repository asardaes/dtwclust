#include "distmat.h"

#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>

#include "concrete-distmat.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<Distmat>
DistmatFactory::create(const SEXP& MAT_TYPE, const SEXP& D)
{
    std::string type = Rcpp::as<std::string>(MAT_TYPE);
    if (type == "R_MATRIX") {
        return std::make_shared<RDistmat>(D);
    }
    else {
        Rcpp::stop("Unknown matrix type"); // nocov
    }
}

} // namespace dtwclust
