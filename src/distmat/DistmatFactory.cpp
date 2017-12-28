#include <RcppArmadillo.h>
#include <memory> // make_shared
#include <string>

#include "distmat.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<Distmat>
DistmatFactory::create(const SEXP& MAT_TYPE, const SEXP& D)
{
    string type = Rcpp::as<string>(MAT_TYPE);
    if (type == "R_MATRIX") {
        return std::make_shared<RDistmat>(D);
    }
    else if (type == "BIG_MATRIX") {
        return std::make_shared<BigmemoryDistmat>(D);
    }
    else {
        Rcpp::stop("Unknown matrix type"); // nocov
    }
}

} // namespace dtwclust
