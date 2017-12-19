#include <Rcpp.h>
#include <algorithm> // std::max
#include <math.h> // exp, log, pow
#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* point-wise squared Euclidean norm */
// =================================================================================================

// univariate
double inline squared_euclidean(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
                                const int i, const int j)
{
    return pow(x[i] - y[j], 2);
}

// multivariate
double inline squared_euclidean(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y,
                                const int i, const int j)
{
    int num_col = x.ncol(); // must be same as y.ncol(), checked in R
    double d = 0;
    for (int k = 0; k < num_col; k++) d += pow(x(i,k) - y(j,k), 2);
    return d;
}

// =================================================================================================
/* soft min operator */
// =================================================================================================

double inline soft_min(double a, double b, double c, const double gamma)
{
    a /= -gamma;
    b /= -gamma;
    c /= -gamma;
    double max_val = std::max(std::max(a, b), c);
    double temp = 0;
    temp += exp(a - max_val);
    temp += exp(b - max_val);
    temp += exp(c - max_val);
    return -gamma * (log(temp) + max_val);
}

// =================================================================================================
/* dynamic programming recursion */
// =================================================================================================

// univariate
double dp_recursion(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
                    Rcpp::NumericMatrix& costmat, const double gamma,
                    SEXP& DISTMAT)
{
    bool save_norm = !Rf_isNull(DISTMAT);
    for (int i = 1; i <= x.length(); i++) {
        for (int j = 1; j <= y.length(); j++)
        {
            double point_norm = squared_euclidean(x, y, i-1, j-1);
            costmat(i,j) = point_norm + soft_min(costmat(i-1, j),
                                                 costmat(i-1, j-1),
                                                 costmat(i, j-1),
                                                 gamma);
            if (save_norm) {
                Rcpp::NumericMatrix dm(DISTMAT);
                dm(i-1, j-1) = point_norm;
            }
        }
    }
    return costmat(x.length(), y.length());
}

// multivariate
double dp_recursion(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y,
                    Rcpp::NumericMatrix& costmat, const double gamma,
                    SEXP& DISTMAT)
{
    bool save_norm = !Rf_isNull(DISTMAT);
    for (int i = 1; i <= x.nrow(); i++) {
        for (int j = 1; j <= y.nrow(); j++)
        {
            double point_norm = squared_euclidean(x, y, i-1, j-1);
            costmat(i,j) = point_norm + soft_min(costmat(i-1, j),
                                                 costmat(i-1, j-1),
                                                 costmat(i, j-1),
                                                 gamma);
            if (save_norm) {
                Rcpp::NumericMatrix dm(DISTMAT);
                dm(i-1, j-1) = point_norm;
            }
        }
    }
    return costmat(x.nrow(), y.nrow());
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP DISTMAT, SEXP MV) {
    BEGIN_RCPP
    bool is_multivariate = Rcpp::as<bool>(MV);
    double gamma = Rcpp::as<double>(GAMMA);
    Rcpp::NumericMatrix costmat(COSTMAT);
    // initialize costmat values
    costmat(0,0) = 0;
    for (int i = 1; i < costmat.nrow(); i++) costmat(i,0) = R_PosInf;
    for (int j = 1; j < costmat.ncol(); j++) costmat(0,j) = R_PosInf;
    // compute distance
    if (is_multivariate) {
        Rcpp::NumericMatrix x(X), y(Y);
        return Rcpp::wrap(dp_recursion(x, y, costmat, gamma, DISTMAT));

    } else {
        Rcpp::NumericVector x(X), y(Y);
        return Rcpp::wrap(dp_recursion(x, y, costmat, gamma, DISTMAT));
    }
    END_RCPP
}

} // namespace dtwclust
