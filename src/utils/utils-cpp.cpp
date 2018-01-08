#include "utils++.h"

#include <RcppArmadillo.h>

namespace dtwclust {

// =================================================================================================
/* Force symmetry helper */
// =================================================================================================

RcppExport SEXP force_lb_symmetry(SEXP X) {
    BEGIN_RCPP
    Rcpp::NumericMatrix matrix(X);
    for (int i = 1; i < matrix.nrow(); i++) {
        Rcpp::checkUserInterrupt();
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

// =================================================================================================
/* helper kahan_sum */
// =================================================================================================

double kahan_sum(const Rcpp::NumericVector& x) {
    double sum = 0, c = 0;
    for (const double& i : x) {
        double y = i - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

} // namespace dtwclust
