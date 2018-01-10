#include "distances++.h"

#include <algorithm> // std::max
#include <math.h> // exp, log, pow
#include <utility> // std::move

#include <RcppArmadillo.h>

#include "../utils/utils.h" // d2s

namespace dtwclust {

// =================================================================================================
/* point-wise squared Euclidean norm */
// =================================================================================================

double inline squared_euclidean(const double * const x, const double * const y,
                                const int i, const int j,
                                const int x_nrows, const int y_nrows, const int ncols)
{
    double d = 0;
    for (int k = 0; k < ncols; k++)
        d += pow(x[d2s(i, k, x_nrows)] - y[d2s(j, k, y_nrows)], 2);
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

template<typename SeriesType, typename MatrixType>
double dp_recursion(const SeriesType& x, const SeriesType&y,
                    MatrixType& costmat, const double gamma,
                    const int nx, const int ny, const int num_vars,
                    const bool save_norm, MatrixType& distmat)
{
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            double point_norm = squared_euclidean(
                &x[0], &y[0], i-1, j-1, nx, ny, num_vars);
            costmat(i,j) = point_norm + soft_min(costmat(i-1, j),
                                                 costmat(i-1, j-1),
                                                 costmat(i, j-1),
                                                 gamma);
            if (save_norm)
                distmat(i-1, j-1) = point_norm;
        }
    }
    return costmat(nx, ny);
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP DISTMAT, SEXP MV) {
    BEGIN_RCPP
    Rcpp::NumericMatrix costmat(COSTMAT), distmat;
    bool is_multivariate = Rcpp::as<bool>(MV);
    double gamma = Rcpp::as<double>(GAMMA);
    bool save_norm = !Rf_isNull(DISTMAT);
    if (save_norm) distmat = std::move(Rcpp::NumericMatrix(DISTMAT));
    // initialize costmat values
    costmat(0,0) = 0;
    for (int i = 1; i < costmat.nrow(); i++) costmat(i,0) = R_PosInf;
    for (int j = 1; j < costmat.ncol(); j++) costmat(0,j) = R_PosInf;
    // compute distance
    if (is_multivariate) {
        Rcpp::NumericMatrix x(X), y(Y);
        return Rcpp::wrap(
            dp_recursion<Rcpp::NumericMatrix, Rcpp::NumericMatrix>(
                    x, y, costmat, gamma, x.nrow(), y.nrow(), x.ncol(), save_norm, distmat));

    } else {
        Rcpp::NumericVector x(X), y(Y);
        return Rcpp::wrap(
            dp_recursion<Rcpp::NumericVector, Rcpp::NumericMatrix>(
                    x, y, costmat, gamma, x.length(), y.length(), 1, save_norm, distmat));
    }
    END_RCPP
}

} // namespace dtwclust
