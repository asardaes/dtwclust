#include <Rcpp.h>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* Force symmetry of square matrix */
// =================================================================================================

RcppExport SEXP force_symmetry(SEXP X, SEXP LOWER) {
    BEGIN_RCPP
    Rcpp::NumericMatrix matrix(X);
    bool lower = Rcpp::as<bool>(LOWER);

    if (lower) {
        for (int i = 1; i < matrix.nrow(); i++) {
            R_CheckUserInterrupt();
            for (int j = 0; j < i; j++) matrix(j,i) = matrix(i,j);
        }
    } else {
        for (int j = 1; j < matrix.ncol(); j++) {
            R_CheckUserInterrupt();
            for (int i = 0; i < j; i++) matrix(j,i) = matrix(i,j);
        }
    }

    return R_NilValue;
    END_RCPP
}

// =================================================================================================
/* Force symmetry of matrix with DTW lower bounds */
// =================================================================================================

RcppExport SEXP force_lb_symmetry(SEXP X) {
    BEGIN_RCPP
    Rcpp::NumericMatrix matrix(X);
    for (int i = 1; i < matrix.nrow(); i++) {
        R_CheckUserInterrupt();
        for (int j = 0; j < i; j++) {
            double lb1 = matrix(i,j);
            double lb2 = matrix(j,i);
            if (lb1 > lb2)
                matrix(j,i) = lb1;
            else
                matrix(i,j) = lb2;
        }
    }
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
