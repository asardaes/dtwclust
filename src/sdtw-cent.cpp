#include <Rcpp.h>
#include <math.h> // exp
#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* backward recursions
 *   Includes the calculation of soft-DTW gradient + jacobian product.
 */
// =================================================================================================

// univariate
void update_gradient(Rcpp::NumericVector& gradient, const double weight, const double gamma,
                     Rcpp::NumericMatrix& cm, Rcpp::NumericMatrix& dm, Rcpp::NumericMatrix& em,
                     const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
    int m = x.length(), n = y.length();
    double a, b, c, grad;
    for (int i = 1; i <= m; i++) {
        dm(i-1, n) = 0;
        cm(i, n+1) = R_NegInf;
    }
    for (int j = 1; j <= n; j++) {
        dm(m, j-1) = 0;
        cm(m+1, j) = R_NegInf;
    }
    cm(m+1, n+1) = cm(m,n);
    dm(m,n) = 0;
    em.fill(0);
    em((m+1)%2, n+1) = 1;
    for (int i = m; i > 0; i--) {
        for (int j = n; j > 0; j--) {
            a = exp((cm(i+1, j) - cm(i,j) - dm(i, j-1)) / gamma);
            b = exp((cm(i, j+1) - cm(i,j) - dm(i-1, j)) / gamma);
            c = exp((cm(i+1, j+1) - cm(i,j) - dm(i,j)) / gamma);
            em(i%2, j) = a * em((i+1)%2, j) + b * em(i%2, j+1) + c * em((i+1)%2, j+1);
        }
        grad = 0;
        for (int jj = 0; jj < n; jj++) grad += em(i%2, jj+1) * 2 * (x[i-1] - y[jj]);
        gradient[i-1] += weight * grad;
        if (i == m) em((m+1)%2, n+1) = 0;
    }
}

// multivariate
void update_gradient(Rcpp::NumericMatrix& gradient, const double weight, const double gamma,
                     Rcpp::NumericMatrix& cm, Rcpp::NumericMatrix& dm, Rcpp::NumericMatrix& em,
                     const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y)
{
    int m = x.nrow(), n = y.nrow(), dim = x.ncol();
    Rcpp::NumericVector grad = Rcpp::NumericVector(dim);
    double a, b, c;
    for (int i = 1; i <= m; i++) {
        dm(i-1, n) = 0;
        cm(i, n+1) = R_NegInf;
    }
    for (int j = 1; j <= n; j++) {
        dm(m, j-1) = 0;
        cm(m+1, j) = R_NegInf;
    }
    cm(m+1, n+1) = cm(m,n);
    dm(m,n) = 0;
    em.fill(0);
    em((m+1)%2, n+1) = 1;
    for (int i = m; i > 0; i--) {
        for (int j = n; j > 0; j--) {
            a = exp((cm(i+1, j) - cm(i,j) - dm(i, j-1)) / gamma);
            b = exp((cm(i, j+1) - cm(i,j) - dm(i-1, j)) / gamma);
            c = exp((cm(i+1, j+1) - cm(i,j) - dm(i,j)) / gamma);
            em(i%2, j) = a * em((i+1)%2, j) + b * em(i%2, j+1) + c * em((i+1)%2, j+1);
        }
        grad.fill(0);
        for (int jj = 0; jj < n; jj++) {
            for (int k = 0; k < dim; k++) {
                grad[k] += em(i%2, jj+1) * 2 * (x(i-1, k) - y(jj, k));
            }
        }
        for (int k = 0; k < dim; k++) gradient(i-1, k) += weight * grad[k];
        if (i == m) em((m+1)%2, n+1) = 0;
    }
}

// =================================================================================================
/* main C++ loops */
// =================================================================================================

// univariate
SEXP sdtw_cent_uv(const Rcpp::List& series, const SEXP& CENT,
                  const Rcpp::NumericVector& weights, const SEXP& GAMMA, const SEXP& MULTIVARIATE,
                  SEXP& COSTMAT, SEXP& DISTMAT, SEXP& EM)
{
    Rcpp::NumericMatrix cm(COSTMAT), dm(DISTMAT), em(EM);
    Rcpp::NumericVector cent(CENT);
    Rcpp::NumericVector gradient = Rcpp::clone(cent);
    gradient.fill(0);
    double objective = 0;
    for (int i = 0; i < series.length(); i++) {
        SEXP X = series[i];
        double d = Rcpp::as<double>(soft_dtw(CENT, X, GAMMA, COSTMAT, DISTMAT, MULTIVARIATE));
        objective += weights[i] * d;
        Rcpp::NumericVector x(X);
        // notice in the call below 'cent' is taken as 'x' and what here is 'x' is 'y' there
        update_gradient(gradient, weights[i], Rcpp::as<double>(GAMMA),
                        cm, dm, em,
                        cent, x);
    }
    Rcpp::List ret = Rcpp::List::create(
        Rcpp::_["objective"] = objective,
        Rcpp::_["gradient"] = gradient
    );
    return ret;
}

// multivariate
SEXP sdtw_cent_mv(const Rcpp::List& series, const SEXP& CENT,
                  const Rcpp::NumericVector& weights, const SEXP& GAMMA, const SEXP& MULTIVARIATE,
                  SEXP& COSTMAT, SEXP& DISTMAT, SEXP& EM)
{
    Rcpp::NumericMatrix cm(COSTMAT), dm(DISTMAT), em(EM);
    Rcpp::NumericMatrix cent(CENT);
    Rcpp::NumericMatrix gradient = Rcpp::clone(cent);
    gradient.fill(0);
    double objective = 0;
    for (int i = 0; i < series.length(); i++) {
        SEXP X = series[i];
        double d = Rcpp::as<double>(soft_dtw(CENT, X, GAMMA, COSTMAT, DISTMAT, MULTIVARIATE));
        objective += weights[i] * d;
        Rcpp::NumericMatrix x(X);
        // notice in the call below 'cent' is taken as 'x' and what here is 'x' is 'y' there
        update_gradient(gradient, weights[i], Rcpp::as<double>(GAMMA),
                        cm, dm, em,
                        cent, x);
    }
    Rcpp::List ret = Rcpp::List::create(
        Rcpp::_["objective"] = objective,
        Rcpp::_["gradient"] = gradient
    );
    return ret;
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID, SEXP GAMMA, SEXP WEIGHTS, SEXP MV,
                          SEXP COSTMAT, SEXP DISTMAT, SEXP EM)
{
    BEGIN_RCPP
    // compute objective and gradient
    if (Rcpp::as<bool>(MV))
        return sdtw_cent_mv(SERIES, CENTROID, WEIGHTS, GAMMA, MV, COSTMAT, DISTMAT, EM);
    else
        return sdtw_cent_uv(SERIES, CENTROID, WEIGHTS, GAMMA, MV, COSTMAT, DISTMAT, EM);
    END_RCPP
}

} // namespace dtwclust
