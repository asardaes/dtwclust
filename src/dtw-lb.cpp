#include <Rcpp.h>
#include "dtwclustpp.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* dtw_basic proxy */
// =================================================================================================

double dtwb(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& dots)
{
    SEXP X = PROTECT(Rcpp::wrap(x));
    SEXP Y = PROTECT(Rcpp::wrap(y));
    SEXP NX = PROTECT(Rcpp::wrap(x.length()));
    SEXP NY = PROTECT(Rcpp::wrap(y.length()));
    SEXP NV = PROTECT(Rcpp::wrap(1));

    SEXP window = dots["window.size"];
    SEXP norm = dots["norm"];
    SEXP step = dots["step.pattern"];
    SEXP backtrack = dots["backtrack"];
    SEXP gcm = dots["gcm"];

    double d = Rcpp::as<double>(dtw_basic(X, Y, window, NX, NY, NV, norm, step, backtrack, gcm));
    UNPROTECT(5);
    return d;
}

// =================================================================================================
/* find nearest neighbors row-wise */
// =================================================================================================

void set_nn(const Rcpp::NumericMatrix& distmat, Rcpp::IntegerVector& nn)
{
    for (int i = 0; i < distmat.nrow(); i++) {
        double d = distmat(i, 0);
        nn[i] = 0;
        for (int j = 1; j < distmat.ncol(); j++) {
            double temp = distmat(i, j);
            if (temp < d) {
                d = temp;
                nn[i] = j;
            }
        }
    }
}

// =================================================================================================
/* check if updates are finished based on indices */
// =================================================================================================

bool check_finished(const Rcpp::IntegerVector& nn, const Rcpp::IntegerVector& nn_prev,
                    Rcpp::LogicalVector& changed)
{
    bool finished = true;
    for (int i = 0; i < nn.length(); i++) {
        if (nn[i] != nn_prev[i]) {
            changed[i] = true;
            finished = false;

        } else {
            changed[i] = false;
        }
    }
    return finished;
}

// =================================================================================================
/* main C++ function */
// =================================================================================================

void dtw_lb_cpp(const Rcpp::List& x, const Rcpp::List& y,
                Rcpp::NumericMatrix& distmat,
                const Rcpp::List& dots)
{
    Rcpp::IntegerVector id_nn(distmat.nrow()), id_nn_prev(distmat.nrow());
    Rcpp::LogicalVector id_changed(distmat.nrow());

    set_nn(distmat, id_nn);
    for (int i = 0; i < id_nn.length(); i++) id_nn_prev[i] = id_nn[i] + 1; // initialize different

    while (!check_finished(id_nn, id_nn_prev, id_changed)) {
        // update nn_prev
        for (int i = 0; i < id_nn.length(); i++) id_nn_prev[i] = id_nn[i];

        // calculate dtw distance if necessary
        for (int i = 0; i < id_changed.length(); i++) {
            if (id_changed[i]) {
                int j = id_nn[i];
                Rcpp::NumericVector this_x = x(i), this_y = y(j);
                distmat(i, j) = dtwb(this_x, this_y, dots);
            }
        }

        // update nearest neighbors
        set_nn(distmat, id_nn);
    }
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP DOTS)
{
BEGIN_RCPP
    Rcpp::List x(X), y(Y), dots(DOTS);
    Rcpp::NumericMatrix distmat(D);
    dtw_lb_cpp(x, y, distmat, dots);
    return R_NilValue;
END_RCPP
}

} // namespace dtwclust
