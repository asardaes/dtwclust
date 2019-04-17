// for whatever reason, OSX doesn't like if I include distmat.h first
// probably due to RcppParallel appearing before Rcpp?

#include <memory> // make_shared
#include <string>

#include <Rcpp.h>
#include <RcppParallel.h>

#include "distmat.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<Distmat> DistmatFactory::create(const SEXP& MAT_TYPE, const SEXP& D) {
    std::string type = Rcpp::as<std::string>(MAT_TYPE);
    if (type == "R_MATRIX") {
        return std::make_shared<RDistmat>(D);
    }
    else {
        Rcpp::stop("Unknown matrix type"); // nocov
    }
}

// =================================================================================================
/* R matrix distmat */
// =================================================================================================

// constructor
RDistmat::RDistmat(const SEXP& D)
    : distmat_(RcppParallel::RMatrix<double>(Rcpp::NumericMatrix(D)))
{ }

// operator() for assignment
double& RDistmat::operator() (const id_t i, const id_t j) {
    return distmat_(i,j);
}

// dimensions
id_t RDistmat::nrow() const { return distmat_.nrow(); }
id_t RDistmat::ncol() const { return distmat_.ncol(); }

} // namespace dtwclust
