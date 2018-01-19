#include "distances.h"

#include <algorithm> // std::max
#include <math.h> // exp, log, pow

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

// template used by the gateway function
template<typename SeriesType>
double dp_recursion(const SeriesType& x, const SeriesType& y,
                    Rcpp::NumericMatrix& costmat, const double gamma,
                    const int nx, const int ny, const int num_vars)
{
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            double point_norm = squared_euclidean(&x[0], &y[0], i-1, j-1, nx, ny, num_vars);
            costmat(i,j) = point_norm + soft_min(costmat(i-1, j),
                                                 costmat(i-1, j-1),
                                                 costmat(i, j-1),
                                                 gamma);
        }
    }
    return costmat(nx, ny);
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV)
{
    BEGIN_RCPP
    Rcpp::NumericMatrix costmat(COSTMAT);
    bool is_multivariate = Rcpp::as<bool>(MV);
    double gamma = Rcpp::as<double>(GAMMA);
    // initialize costmat values
    costmat(0,0) = 0;
    for (int i = 1; i < costmat.nrow(); i++) costmat(i,0) = R_PosInf;
    for (int j = 1; j < costmat.ncol(); j++) costmat(0,j) = R_PosInf;
    // compute distance
    if (is_multivariate) {
        Rcpp::NumericMatrix x(X), y(Y);
        return Rcpp::wrap(
            dp_recursion<Rcpp::NumericMatrix>(
                    x, y, costmat, gamma, x.nrow(), y.nrow(), x.ncol()));
    }
    else {
        Rcpp::NumericVector x(X), y(Y);
        return Rcpp::wrap(
            dp_recursion<Rcpp::NumericVector>(
                    x, y, costmat, gamma, x.length(), y.length(), 1));
    }
    END_RCPP
}

// =================================================================================================
/* thread-safe version for the distance calculator */
// =================================================================================================

double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, double * const costmat,
            const bool save_norm, double * const distmat)
{
    // initialize costmat values
    costmat[0] = 0;
    for (int i = 1; i < nx+2; i++) costmat[d2s(i,0,nx+2)] = R_PosInf;
    for (int j = 1; j < ny+2; j++) costmat[d2s(0,j,nx+2)] = R_PosInf;
    // compute distance
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            double point_norm = squared_euclidean(x, y, i-1, j-1, nx, ny, num_vars);
            costmat[d2s(i,j,nx+2)] = point_norm + soft_min(costmat[d2s(i-1, j, nx+2)],
                                                           costmat[d2s(i-1, j-1, nx+2)],
                                                           costmat[d2s(i, j-1, nx+2)],
                                                           gamma);
            if (save_norm)
                distmat[d2s(i-1,j-1,nx+1)] = point_norm;
        }
    }
    return costmat[d2s(nx,ny,nx+2)];
}

} // namespace dtwclust
