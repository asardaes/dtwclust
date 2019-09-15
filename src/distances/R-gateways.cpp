#include "R-gateways.h"

#include <Rcpp.h>

#include "details.h"
#include "../utils/SurrogateMatrix.h"
#include "../utils/utils.h" // id_t

namespace dtwclust {

// =================================================================================================
/* dtw_basic */
// =================================================================================================

extern "C" SEXP dtw_basic(SEXP X, SEXP Y, SEXP WINDOW,
                          SEXP X_LEN, SEXP Y_LEN, SEXP NUM_VAR,
                          SEXP NORM, SEXP STEP, SEXP BACKTRACK, SEXP NORMALIZE, SEXP SQRT_DIST,
                          SEXP LCM)
{
    BEGIN_RCPP
    double *lcm = REAL(LCM);
    int nx = Rf_asInteger(X_LEN);
    int ny = Rf_asInteger(Y_LEN);
    int nv = Rf_asInteger(NUM_VAR);
    SurrogateMatrix<const double> x(nx, nv, REAL(X));
    SurrogateMatrix<const double> y(ny, nv, REAL(Y));

    if (Rf_asLogical(BACKTRACK)) {
        SurrogateMatrix<double> wrapped_lcm(nx + 1, ny + 1, lcm);

        // longest possible path, length will be adjusted in R
        SEXP INDEX1 = PROTECT(Rf_allocVector(INTSXP, nx + ny));
        SEXP INDEX2 = PROTECT(Rf_allocVector(INTSXP, nx + ny));
        SurrogateMatrix<int> index1(nx + ny, 1, INTEGER(INDEX1));
        SurrogateMatrix<int> index2(nx + ny, 1, INTEGER(INDEX2));
        int path = 0;

        // calculate distance
        double d = dtw_basic(wrapped_lcm, x, y,
                             Rf_asInteger(WINDOW), Rf_asReal(NORM), Rf_asReal(STEP),
                             Rf_asLogical(NORMALIZE), Rf_asLogical(SQRT_DIST),
                             index1, index2, path);

        // put results in a list
        SEXP list_names = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(list_names, 0, Rf_mkChar("distance"));
        SET_STRING_ELT(list_names, 1, Rf_mkChar("index1"));
        SET_STRING_ELT(list_names, 2, Rf_mkChar("index2"));
        SET_STRING_ELT(list_names, 3, Rf_mkChar("path"));

        SEXP ret = PROTECT(Rf_allocVector(VECSXP, 4));
        SET_VECTOR_ELT(ret, 0, PROTECT(Rf_ScalarReal(d)));
        SET_VECTOR_ELT(ret, 1, INDEX1);
        SET_VECTOR_ELT(ret, 2, INDEX2);
        SET_VECTOR_ELT(ret, 3, PROTECT(Rf_ScalarInteger(path)));
        Rf_setAttrib(ret, R_NamesSymbol, list_names);

        UNPROTECT(6);
        return ret;
    }
    else {
        SurrogateMatrix<double> wrapped_lcm(2, ny + 1, lcm);

        // calculate distance
        double d = dtw_basic(wrapped_lcm, x, y,
                             Rf_asInteger(WINDOW), Rf_asReal(NORM), Rf_asReal(STEP),
                             Rf_asLogical(NORMALIZE), Rf_asLogical(SQRT_DIST));

        return Rcpp::wrap(d);
    }
    END_RCPP
}

// =================================================================================================
/* lbi */
// =================================================================================================

extern "C" SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U)
{
    BEGIN_RCPP
    Rcpp::NumericVector x(X), y(Y), lower(L), upper(U);
    SurrogateMatrix<const double> temp_x(x.length(), 1, &x[0]);
    SurrogateMatrix<const double> temp_y(y.length(), 1, &y[0]);
    SurrogateMatrix<const double> temp_l(lower.length(), 1, &lower[0]);
    SurrogateMatrix<const double> temp_u(upper.length(), 1, &upper[0]);
    SurrogateMatrix<double> L2(x.length(), 1);
    SurrogateMatrix<double> U2(x.length(), 1);
    SurrogateMatrix<double> H(x.length(), 1);
    SurrogateMatrix<double> LB(x.length(), 1);
    return Rcpp::wrap(lbi_core(temp_x, temp_y,
                               Rcpp::as<unsigned int>(WINDOW), Rcpp::as<int>(P),
                               temp_l, temp_u, L2, U2, H, LB));
    END_RCPP
}

// =================================================================================================
/* lbk */
// =================================================================================================

extern "C" SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U)
{
    BEGIN_RCPP
    Rcpp::NumericVector x(X), lower(L), upper(U);
    SurrogateMatrix<const double> temp_x(x.length(), 1, &x[0]);
    SurrogateMatrix<const double> temp_l(lower.length(), 1, &lower[0]);
    SurrogateMatrix<const double> temp_u(upper.length(), 1, &upper[0]);
    SurrogateMatrix<double> H(x.length(), 1);
    return Rcpp::wrap(lbk_core(temp_x, Rcpp::as<int>(P), temp_l, temp_u, H));
    END_RCPP
}

// =================================================================================================
/* logGAK */
// =================================================================================================

extern "C" SEXP logGAK(SEXP X, SEXP Y, SEXP NX, SEXP NY, SEXP NUM_VAR,
                       SEXP SIGMA, SEXP WINDOW, SEXP LOGS)
{
    BEGIN_RCPP
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

    int nx = Rf_asInteger(NX);
    int ny = Rf_asInteger(NY);
    int num_var = Rf_asInteger(NUM_VAR);
    int nlogs = (nx > ny) ? nx + 1 : ny + 1;
    SurrogateMatrix<const double> x(nx, num_var, REAL(X));
    SurrogateMatrix<const double> y(ny, num_var, REAL(Y));
    SurrogateMatrix<double> logs(nlogs, 3, REAL(LOGS));
    return Rf_ScalarReal(
        logGAK_c(x, y, Rf_asReal(SIGMA), static_cast<id_t>(Rf_asInteger(WINDOW)), logs)
    );
    END_RCPP
}

// =================================================================================================
/* soft-DTW */
// =================================================================================================

extern "C" SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV)
{
    BEGIN_RCPP
    Rcpp::NumericMatrix costmat(COSTMAT);
    SurrogateMatrix<double> cm(costmat.nrow(), costmat.ncol(), &costmat[0]);
    bool is_multivariate = Rcpp::as<bool>(MV);
    double gamma = Rcpp::as<double>(GAMMA);

    // compute distance
    if (is_multivariate) {
        Rcpp::NumericMatrix x(X), y(Y);
        SurrogateMatrix<const double> temp_x(x.nrow(), x.ncol(), &x[0]);
        SurrogateMatrix<const double> temp_y(y.nrow(), y.ncol(), &y[0]);
        return Rcpp::wrap(sdtw(temp_x, temp_y, gamma, cm));
    }
    else {
        Rcpp::NumericVector x(X), y(Y);
        SurrogateMatrix<const double> temp_x(x.length(), 1, &x[0]);
        SurrogateMatrix<const double> temp_y(y.length(), 1, &y[0]);
        return Rcpp::wrap(sdtw(temp_x, temp_y, gamma, cm));
    }
    END_RCPP
}

} // namespace dtwclust
