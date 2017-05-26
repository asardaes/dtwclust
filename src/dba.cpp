#include <Rcpp.h>
#include <cmath>
#include "dtwclustpp.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* shared variables */
// =================================================================================================

SEXP window, norm, step, backtrack;
Rcpp::List series, dots;
Rcpp::IntegerVector index1, index2;
Rcpp::NumericMatrix gcm;
int max_iter, nx, ny, nv, begin;
double delta;
bool trace;

// =================================================================================================
/* set alignment with dtw_basic (yes, it looks ugly) */
// =================================================================================================

void uv_set_alignment(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
    SEXP X = PROTECT(Rcpp::wrap(x));
    SEXP Y = PROTECT(Rcpp::wrap(y));
    SEXP GCM = PROTECT(Rcpp::wrap(gcm));
    SEXP NX = PROTECT(Rcpp::wrap(nx));
    SEXP NY = PROTECT(Rcpp::wrap(ny));
    SEXP NV = PROTECT(Rcpp::wrap(nv));

    Rcpp::List alignment(dtw_basic(X, Y, window, NX, NY, NV, norm, step, backtrack, GCM));

    index1 = alignment["index1"];
    index2 = alignment["index2"];
    begin = alignment["path"];
    UNPROTECT(6);
}

void mv_set_alignment(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y)
{
    SEXP X = PROTECT(Rcpp::wrap(x));
    SEXP Y = PROTECT(Rcpp::wrap(y));
    SEXP GCM = PROTECT(Rcpp::wrap(gcm));
    SEXP NX = PROTECT(Rcpp::wrap(nx));
    SEXP NY = PROTECT(Rcpp::wrap(ny));
    SEXP NV = PROTECT(Rcpp::wrap(nv));

    Rcpp::List alignment(dtw_basic(X, Y, window, NX, NY, NV, norm, step, backtrack, GCM));

    index1 = alignment["index1"];
    index2 = alignment["index2"];
    begin = alignment["path"];
    UNPROTECT(6);
}

// =================================================================================================
/* sum step for vectors and matrices */
// =================================================================================================

void uv_sum_step(const Rcpp::NumericVector& x, Rcpp::NumericVector& cent, Rcpp::IntegerVector& n,
                 Rcpp::NumericVector& kahan_c, Rcpp::NumericVector& kahan_y, Rcpp::NumericVector& kahan_t)
{
    for (int i = begin - 1; i >= 0; i--) {
        int i1 = index1[i] - 1;
        int i2 = index2[i] - 1;
        kahan_y[i2] = x[i1] - kahan_c[i2];
        kahan_t[i2] = cent[i2] + kahan_y[i2];
        kahan_c[i2] = (kahan_t[i2] - cent[i2]) - kahan_y[i2];
        cent[i2] = kahan_t[i2];
        n[i2] += 1;
    }
}

void mv_sum_step(const Rcpp::NumericMatrix& x, Rcpp::NumericMatrix& cent, Rcpp::IntegerMatrix& n,
                 Rcpp::NumericMatrix& kahan_c, Rcpp::NumericMatrix& kahan_y, Rcpp::NumericMatrix& kahan_t,
                 int col)
{
    for (int i = begin - 1; i >= 0; i--) {
        int i1 = index1[i] - 1;
        int i2 = index2[i] - 1;
        kahan_y(i2, col) = x(i1, col) - kahan_c(i2, col);
        kahan_t(i2, col) = cent(i2, col) + kahan_y(i2, col);
        kahan_c(i2, col) = (kahan_t(i2, col) - cent(i2, col)) - kahan_y(i2, col);
        cent(i2, col) = kahan_t(i2, col);
        n(i2, col) += 1;
    }
}

// =================================================================================================
/* average step for vectors and matrices */
// =================================================================================================

bool uv_average_step(Rcpp::NumericVector& new_cent,
                     Rcpp::IntegerVector& num_vals,
                     Rcpp::NumericVector& ref_cent)
{
    bool converged = true;
    for (int i = 0; i < ny; i++) {
        new_cent[i] /= num_vals[i];
        if (std::abs(new_cent[i] - ref_cent[i]) >= delta) converged = false;
        ref_cent[i] = new_cent[i];
    }

    return converged;
}

bool mv_average_step(Rcpp::NumericMatrix& new_cent,
                     Rcpp::IntegerMatrix& num_vals,
                     Rcpp::NumericMatrix& ref_cent)
{
    bool converged = true;
    for (int j = 0; j < nv; j++) {
        for (int i = 0; i < ny; i++) {
            new_cent(i, j) /= num_vals(i, j);
            if (std::abs(new_cent(i, j) - ref_cent(i, j)) >= delta) converged = false;
            ref_cent(i, j) = new_cent(i, j);
        }
    }

    return converged;
}

// =================================================================================================
/* helper functions for all */
// =================================================================================================

int max_lengths(bool mv)
{
    int max_length = 0;
    for (int i = 0; i < series.length(); i++) {
        int temp;
        if (mv) {
            Rcpp::NumericMatrix x = series[i];
            temp = x.nrow();

        } else {
            Rcpp::NumericVector x = series[i];
            temp = x.length();
        }

        if (temp > max_length) max_length = temp;
    }
    return max_length;
}

void reset_vectors(Rcpp::NumericVector& new_cent,
                   Rcpp::IntegerVector& num_vals,
                   Rcpp::NumericVector& kahan_c)
{
    new_cent.fill(0);
    num_vals.fill(0);
    kahan_c.fill(0);
}

void reset_matrices(Rcpp::NumericMatrix& mat_cent,
                    Rcpp::IntegerMatrix& mat_vals,
                    Rcpp::NumericMatrix& kahan_c)
{
    mat_cent.fill(0);
    mat_vals.fill(0);
    kahan_c.fill(0);
}

void Rflush()
{
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
}

void print_trace(bool converged, int iter)
{
    if (trace) {
        if (converged) {
            Rcpp::Rcout << " " << iter << " - Converged!" << std::endl;

        } else {
            Rcpp::Rcout << " " << iter << ",";
            if (iter % 10 == 0) Rcpp::Rcout << "\n\t\t";
            Rflush();
        }
    }
}

// =================================================================================================
/* univariate DBA */
// =================================================================================================

SEXP dba_uv(const SEXP& centroid)
{
    Rcpp::NumericVector ref_cent = Rcpp::clone(Rcpp::NumericVector(centroid));
    ny = ref_cent.length();
    nv = 1;

    Rcpp::NumericVector new_cent(ny);
    Rcpp::IntegerVector num_vals(ny);

    // for Kahan summation
    Rcpp::NumericVector kahan_c(ny), kahan_y(ny), kahan_t(ny);

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        reset_vectors(new_cent, num_vals, kahan_c);

        // sum step
        for (int i = 0; i < series.length(); i++) {
            Rcpp::NumericVector x = series[i];
            nx = x.length();
            uv_set_alignment(x, ref_cent);
            uv_sum_step(x, new_cent, num_vals, kahan_c, kahan_y, kahan_t);
        }

        // average step with check for 'convergence' and update ref_cent
        bool converged = uv_average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
    }

    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'" << std::endl;
        Rflush();
    }

    UNPROTECT(4);
    return Rcpp::wrap(new_cent);
}

// =================================================================================================
/* multivariate DBA considering each variable separately */
// =================================================================================================

SEXP dba_mv_by_variable(const SEXP& centroid)
{
    Rcpp::NumericMatrix mv_ref_cent(centroid);
    ny = mv_ref_cent.nrow();
    nv = 1; // careful! this is used by the uv_* functions so it must be 1

    Rcpp::NumericMatrix mat_cent(ny, mv_ref_cent.ncol());
    Rcpp::NumericVector x(max_lengths(true));
    Rcpp::NumericVector ref_cent(ny);
    Rcpp::NumericVector new_cent(ny);
    Rcpp::IntegerVector num_vals(ny);

    // for Kahan summation
    Rcpp::NumericVector kahan_c(ny), kahan_y(ny), kahan_t(ny);

    for (int j = 0; j < mv_ref_cent.ncol(); j++) {
        if (trace) Rcpp::Rcout << "\tDBA Iteration:";
        for (int k = 0; k < ny; k++) ref_cent[k] = mv_ref_cent(k, j);
        int iter = 1;
        while (iter <= max_iter) {
            reset_vectors(new_cent, num_vals, kahan_c);

            // sum step
            for (int i = 0; i < series.length(); i++) {
                Rcpp::NumericMatrix mv_x = series[i];
                nx = mv_x.nrow();
                for (int k = 0; k < nx; k++) x[k] = mv_x(k, j);
                uv_set_alignment(x, ref_cent);
                uv_sum_step(x, new_cent, num_vals, kahan_c, kahan_y, kahan_t);
            }

            // average step with check for 'convergence' and update ref_cent
            bool converged = uv_average_step(new_cent, num_vals, ref_cent);
            print_trace(converged, iter);
            if (converged) break;
            iter++;
        }

        if (iter > max_iter && trace) {
            Rcpp::Rcout << " Did not 'converge'" << std::endl;
            Rflush();
        }

        mat_cent(Rcpp::_, j) = new_cent;
    }

    UNPROTECT(4);
    return Rcpp::wrap(mat_cent);
}

// =================================================================================================
/* multivariate DBA considering each series as a whole */
// =================================================================================================

SEXP dba_mv_by_series(const SEXP& centroid)
{
    Rcpp::NumericMatrix ref_cent = Rcpp::clone(Rcpp::NumericMatrix(centroid));
    ny = ref_cent.nrow();
    nv = ref_cent.ncol();

    Rcpp::NumericMatrix mat_cent(ny, nv);
    Rcpp::IntegerMatrix mat_vals(ny, nv);

    // for Kahan summation
    Rcpp::NumericMatrix kahan_c(ny, nv), kahan_y(ny, nv), kahan_t(ny, nv);

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        reset_matrices(mat_cent, mat_vals, kahan_c);

        // sum step
        for (int i = 0; i < series.length(); i++) {
            Rcpp::NumericMatrix x = series[i];
            nx = x.nrow();
            mv_set_alignment(x, ref_cent);
            for (int j = 0; j < nv; j++)
                mv_sum_step(x, mat_cent, mat_vals, kahan_c, kahan_y, kahan_t, j);
        }

        // average step with check for 'convergence' and update ref_cent
        bool converged = mv_average_step(mat_cent, mat_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
    }

    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'" << std::endl;
        Rflush();
    }

    UNPROTECT(4);
    return Rcpp::wrap(mat_cent);
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP dba(SEXP X, SEXP centroid,
                    SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
                    SEXP multivariate, SEXP mv_ver, SEXP DOTS)
{
BEGIN_RCPP
    series = Rcpp::List(X);

    max_iter = Rcpp::as<int>(MAX_ITER);
    delta = Rcpp::as<double>(DELTA);
    trace = Rcpp::as<bool>(TRACE);

    dots = Rcpp::List(DOTS);
    window = PROTECT(Rcpp::wrap(dots["window.size"]));
    norm = PROTECT(Rcpp::wrap(dots["norm"]));
    step = PROTECT(Rcpp::wrap(dots["step.pattern"]));
    backtrack = PROTECT(Rcpp::wrap(true));

    // doing it in one step complains about ambiguous operator=
    Rcpp::NumericMatrix GCM = dots["gcm"];
    gcm = GCM;

    if (Rcpp::as<bool>(multivariate)) {
        if (Rcpp::as<int>(mv_ver) == 1)
            return dba_mv_by_variable(centroid);
        else
            return dba_mv_by_series(centroid);

    } else return dba_uv(centroid);
END_RCPP
}

} // namespace dtwclust
