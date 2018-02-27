#include "distances.h"

#include <Rcpp.h>

#include "distances-details.h"

namespace dtwclust {

// =================================================================================================
/* dtw_basic */
// =================================================================================================

extern "C" SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
                          SEXP m, SEXP n, SEXP num_var,
                          SEXP norm, SEXP step, SEXP backtrack, SEXP normalize,
                          SEXP distmat)
{
    double d;
    int nx = Rf_asInteger(m);
    int ny = Rf_asInteger(n);
    double* D = REAL(distmat);
    dtwclust_tuple_t tuple[3];

    if (Rf_asLogical(backtrack)) {
        // longest possible path, length will be adjusted in R
        SEXP index1 = PROTECT(Rf_allocVector(INTSXP, nx + ny));
        SEXP index2 = PROTECT(Rf_allocVector(INTSXP, nx + ny));

        // calculate distance
        d = dtw_basic_c(D, tuple,
                        REAL(x), REAL(y), Rf_asInteger(window),
                        nx, ny, Rf_asInteger(num_var),
                        Rf_asReal(norm), Rf_asReal(step), 1);
        if (Rf_asLogical(normalize)) d /= nx + ny;

        // actual length of path
        int path = backtrack_steps(D, nx, ny, INTEGER(index1), INTEGER(index2));

        // put results in a list
        SEXP list_names = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(list_names, 0, Rf_mkChar("distance"));
        SET_STRING_ELT(list_names, 1, Rf_mkChar("index1"));
        SET_STRING_ELT(list_names, 2, Rf_mkChar("index2"));
        SET_STRING_ELT(list_names, 3, Rf_mkChar("path"));

        SEXP ret = PROTECT(Rf_allocVector(VECSXP, 4));
        SET_VECTOR_ELT(ret, 0, PROTECT(Rf_ScalarReal(d)));
        SET_VECTOR_ELT(ret, 1, index1);
        SET_VECTOR_ELT(ret, 2, index2);
        SET_VECTOR_ELT(ret, 3, PROTECT(Rf_ScalarInteger(path)));
        Rf_setAttrib(ret, R_NamesSymbol, list_names);

        UNPROTECT(6);
        return ret;
    }
    else {
        // calculate distance
        d = dtw_basic_c(D, tuple,
                        REAL(x), REAL(y), Rf_asInteger(window),
                        nx, ny, Rf_asInteger(num_var),
                        Rf_asReal(norm), Rf_asReal(step), 0);
        if (Rf_asLogical(normalize)) d /= nx + ny;

        SEXP ret = PROTECT(Rf_ScalarReal(d));
        UNPROTECT(1);
        return ret;
    }
}

// =================================================================================================
/* lbi */
// =================================================================================================

extern "C" SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U)
{
    BEGIN_RCPP
    Rcpp::NumericVector x(X), y(Y), lower(L), upper(U);
    Rcpp::NumericVector L2(x.length()), U2(x.length()), H(x.length()), LB(x.length());
    return Rcpp::wrap(lbi_core(&x[0], &y[0], x.length(),
                               Rcpp::as<unsigned int>(WINDOW), Rcpp::as<int>(P),
                               &lower[0], &upper[0], &L2[0], &U2[0], &H[0], &LB[0]));
    END_RCPP
}

// =================================================================================================
/* lbk */
// =================================================================================================

extern "C" SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U)
{
    BEGIN_RCPP
    Rcpp::NumericVector x(X), lower(L), upper(U);
    Rcpp::NumericVector H(x.length());
    return Rcpp::wrap(lbk_core(&x[0], x.length(), Rcpp::as<int>(P), &lower[0], &upper[0], &H[0]));
    END_RCPP
}

// =================================================================================================
/* logGAK */
// =================================================================================================

extern "C" SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var,
                       SEXP sigma, SEXP window, SEXP logs)
{
    /*
     * Inputs are, in this order
     * A N1 x d matrix (d-variate time series with N1 observations)
     * A N2 x d matrix (d-variate time series with N2 observations)
     * The length of series N1 and N2, and the number of variables d
     * A sigma > 0 parameter for the Gaussian kernel's width
     * A triangular integer which parameterizes the band.
     *   - when triangular = 0, the triangular kernel is not used, the evaluation is on the sum of all
     *   paths, and the complexity if of the order of N1 x N2
     *   - when triangular > 0, the triangular kernel is used and downweights some of the paths that
     *   lay far from the diagonal.
     * A (max(N1,N2) + 1) x 3 matrix of doubles
     */

    int triangular = Rf_asInteger(window);
    int nX = Rf_asInteger(nx);
    int nY = Rf_asInteger(ny);
    double d;

    // If triangular is smaller than the difference in length of the time series,
    // the kernel is equal to zero,
    // i.e. its log is set to -Inf
    if (triangular > 0 && abs(nX - nY) > triangular)
        d = R_NegInf;
    else
        d = logGAK_c(REAL(x), REAL(y),
                     nX, nY, Rf_asInteger(num_var),
                     Rf_asReal(sigma), triangular,
                     REAL(logs));
    return Rf_ScalarReal(d);
}

// =================================================================================================
/* soft-DTW */
// =================================================================================================

// dynamic programming
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

// gateway
extern "C" SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV)
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

} // namespace dtwclust
