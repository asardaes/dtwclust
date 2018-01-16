#include "centroids.h"

#include <algorithm> // std::fill
#include <cmath> // std::abs
#include <cstddef> // std::size_t
#include <utility> // std::move

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distance-calculators/TSTSList.h"
#include "../distances/distances.h" // dtw_basic_par
#include "../utils/utils++.h" // get_grain
#include "../utils/utils.h" // Rflush, d2s

namespace dtwclust {

// =================================================================================================
/* a thread-safe DTW calculator that includes backtracking */
// =================================================================================================

class DtwBacktrackCalculator : public DistanceCalculator
{
public:
    // constructor
    DtwBacktrackCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y, const bool mv)
        : is_multivariate_(mv)
    {
        Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
        window_ = Rcpp::as<int>(dist_args["window.size"]);
        norm_ = Rcpp::as<double>(dist_args["norm"]);
        step_ = Rcpp::as<double>(dist_args["step.pattern"]);
        normalize_ = Rcpp::as<bool>(dist_args["normalize"]);
        if (is_multivariate_) {
            x_mv_ = TSTSList<Rcpp::NumericMatrix>(x);
            y_mv_ = TSTSList<Rcpp::NumericMatrix>(y);
        }
        else {
            x_uv_ = TSTSList<Rcpp::NumericVector>(x);
            y_uv_ = TSTSList<Rcpp::NumericVector>(y);
        }
        // set value of max_len_*_
        max_len_x_ = this->maxLength(x, is_multivariate_);
        max_len_y_ = this->maxLength(y, is_multivariate_);
        // make sure pointers are null
        gcm_ = nullptr;
        index1_ = nullptr;
        index2_ = nullptr;
    }

    // destructor
    ~DtwBacktrackCalculator() {
        if (gcm_) delete[] gcm_;
        if (index1_) delete[] index1_;
        if (index2_) delete[] index2_;
    }

    // calculate for given indices (inherited)
    double calculate(const int i, const int j) override {
        if (is_multivariate_)
            return this->calculate(x_mv_[i], y_mv_[j]);
        else
            return this->calculate(x_uv_[i], y_uv_[j]);
    }

    // clone to setup helpers in each thread
    DtwBacktrackCalculator* clone() const override {
        DtwBacktrackCalculator* ptr = new DtwBacktrackCalculator(*this);
        ptr->gcm_ = new double[(max_len_x_ + 1) * (max_len_y_ + 1)];
        ptr->index1_ = new int[max_len_x_ + max_len_y_];
        ptr->index2_ = new int[max_len_x_ + max_len_y_];
        return ptr;
    }

    // limits
    int xLimit() const override {
        return is_multivariate_ ? x_mv_.length() : x_uv_.length();
    }
    int yLimit() const override {
        return is_multivariate_ ? y_mv_.length() : y_uv_.length();
    }

    // input series (univariate)
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    // input series (multivariate)
    TSTSList<Rcpp::NumericMatrix> x_mv_, y_mv_;
    // helpers for backtracking
    int path_, *index1_, *index2_;

private:
    // univariate calculate
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& y)
    {
        if (!gcm_ || !index1_ || !index2_) return -1;
        int nx = x.length();
        int ny = y.length();
        int num_var = 1;
        return dtw_basic_par(&x[0], &y[0],
                             nx, ny, num_var,
                             window_, norm_, step_, normalize_,
                             gcm_, true, index1_, index2_, &path_);
    }

    // multivariate calculate
    double calculate(const RcppParallel::RMatrix<double>& x,
                     const RcppParallel::RMatrix<double>& y)
    {
        if (!gcm_ || !index1_ || !index2_) return -1;
        int nx = x.nrow();
        int ny = y.nrow();
        int num_var = x.ncol();
        return dtw_basic_par(&x[0], &y[0],
                             nx, ny, num_var,
                             window_, norm_, step_, normalize_,
                             gcm_, true, index1_, index2_, &path_);
    }

    // input parameters
    int window_;
    double norm_, step_;
    bool normalize_, is_multivariate_;
    // helper "matrix"
    double* gcm_;
    // to dimension helpers
    int max_len_x_, max_len_y_;
};

// =================================================================================================
/* the parallel worker for the univariate version */
// =================================================================================================

class DbaUv : public RcppParallel::Worker {
public:
    // constructor
    DbaUv(const DtwBacktrackCalculator&& backtrack_calculator,
          const Rcpp::NumericVector& new_cent,
          const Rcpp::IntegerVector& num_vals)
        : backtrack_calculator_(backtrack_calculator)
        , new_cent_(new_cent)
        , num_vals_(num_vals)
        , kahan_c_(new double[new_cent.length()])
        , kahan_y_(new double[new_cent.length()])
        , kahan_t_(new double[new_cent.length()])
    { }

    // destructor
    ~DbaUv() {
        if (kahan_c_) delete[] kahan_c_;
        if (kahan_y_) delete[] kahan_y_;
        if (kahan_t_) delete[] kahan_t_;
    }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();
        // kahan sum step
        std::fill(kahan_c_, kahan_c_ + new_cent_.length(), 0);
        for (std::size_t i = begin; i < end; i++) {
            local_calculator->calculate(i,0);
            RcppParallel::RVector<double>& x = local_calculator->x_uv_[i];
            mutex_.lock();
            for (int ii = local_calculator->path_ - 1; ii >= 0; ii--) {
                int i1 = local_calculator->index1_[ii] - 1;
                int i2 = local_calculator->index2_[ii] - 1;
                kahan_y_[i2] = x[i1] - kahan_c_[i2];
                kahan_t_[i2] = new_cent_[i2] + kahan_y_[i2];
                kahan_c_[i2] = (kahan_t_[i2] - new_cent_[i2]) - kahan_y_[i2];
                new_cent_[i2] = kahan_t_[i2];
                num_vals_[i2] += 1;
            }
            mutex_.unlock();
        }
        // finish
        mutex_.lock();
        delete local_calculator;
        mutex_.unlock();
    }

private:
    const DtwBacktrackCalculator& backtrack_calculator_;
    RcppParallel::RVector<double> new_cent_;
    RcppParallel::RVector<int> num_vals_;
    // sum helpers
    double *kahan_c_, *kahan_y_, *kahan_t_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// =================================================================================================
/* shared variables */
// =================================================================================================

static SEXP window, norm, step, backtrack, normalize, gcm;
static Rcpp::List series;
static Rcpp::IntegerVector index1, index2;
static int max_iter, nx, ny, nv, path, num_threads;
static double delta;
static bool trace;

// =================================================================================================
/* set alignment with dtw_basic (yes, it looks ugly) */
// =================================================================================================

void uv_set_alignment(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y)
{
    SEXP NX = PROTECT(Rcpp::wrap(nx));
    SEXP NY = PROTECT(Rcpp::wrap(ny));
    SEXP NV = PROTECT(Rcpp::wrap(nv));
    Rcpp::List alignment(dtw_basic(x, y, window, NX, NY, NV, norm, step, backtrack, normalize, gcm));
    index1 = alignment["index1"];
    index2 = alignment["index2"];
    path = alignment["path"];
    UNPROTECT(3);
}

void mv_set_alignment(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y)
{
    SEXP NX = PROTECT(Rcpp::wrap(nx));
    SEXP NY = PROTECT(Rcpp::wrap(ny));
    SEXP NV = PROTECT(Rcpp::wrap(nv));
    Rcpp::List alignment(dtw_basic(x, y, window, NX, NY, NV, norm, step, backtrack, normalize, gcm));
    index1 = alignment["index1"];
    index2 = alignment["index2"];
    path = alignment["path"];
    UNPROTECT(3);
}

// =================================================================================================
/* sum step for vectors and matrices */
// =================================================================================================

void uv_sum_step(const Rcpp::NumericVector& x,
                 Rcpp::NumericVector& cent,
                 Rcpp::IntegerVector& n,
                 Rcpp::NumericVector& kahan_c,
                 Rcpp::NumericVector& kahan_y,
                 Rcpp::NumericVector& kahan_t)
{
    for (int i = path - 1; i >= 0; i--) {
        int i1 = index1[i] - 1;
        int i2 = index2[i] - 1;
        kahan_y[i2] = x[i1] - kahan_c[i2];
        kahan_t[i2] = cent[i2] + kahan_y[i2];
        kahan_c[i2] = (kahan_t[i2] - cent[i2]) - kahan_y[i2];
        cent[i2] = kahan_t[i2];
        n[i2] += 1;
    }
}

void mv_sum_step(const Rcpp::NumericMatrix& x,
                 Rcpp::NumericMatrix& cent,
                 Rcpp::IntegerMatrix& n,
                 Rcpp::NumericMatrix& kahan_c,
                 Rcpp::NumericMatrix& kahan_y,
                 Rcpp::NumericMatrix& kahan_t)
{
    for (int j = 0; j < nv; j++) {
        for (int i = path - 1; i >= 0; i--) {
            int i1 = index1[i] - 1;
            int i2 = index2[i] - 1;
            kahan_y(i2, j) = x(i1, j) - kahan_c(i2, j);
            kahan_t(i2, j) = cent(i2, j) + kahan_y(i2, j);
            kahan_c(i2, j) = (kahan_t(i2, j) - cent(i2, j)) - kahan_y(i2, j);
            cent(i2, j) = kahan_t(i2, j);
            n(i2, j) += 1;
        }
    }
}

// =================================================================================================
/* average step with check for 'convergence' and update ref_cent for vectors and matrices */
// =================================================================================================

bool uv_average_step(Rcpp::NumericVector& new_cent,
                     const Rcpp::IntegerVector num_vals,
                     Rcpp::NumericVector& ref_cent)
{
    bool converged = true;
    for (int i = 0; i < ref_cent.length(); i++) {
        new_cent[i] /= num_vals[i];
        if (std::abs(new_cent[i] - ref_cent[i]) >= delta) converged = false;
        ref_cent[i] = new_cent[i];
    }
    return converged;
}

bool mv_average_step(Rcpp::NumericMatrix& new_cent,
                     const Rcpp::IntegerMatrix& num_vals,
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

int max_lengths_mv()
{
    int max_length = 0, temp;
    for (int i = 0; i < series.length(); i++) {
        Rcpp::NumericMatrix x = series[i];
        temp = x.nrow();
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

void print_trace(const bool converged, const int iter)
{
    if (trace) {
        if (converged) {
            Rcpp::Rcout << " " << iter << " - Converged!" << std::endl;

        } else {
            Rcpp::Rcout << " " << iter << ",";
            if (iter % 10 == 0) {
                Rcpp::Rcout << "\n\t\t";
                Rflush();
            }
        }
    }
}

// =================================================================================================
/* univariate DBA */
// =================================================================================================

SEXP dba_uv(const SEXP& X, const Rcpp::NumericVector& centroid,
            const SEXP& DOTS, const bool is_multivariate)
{
    Rcpp::NumericVector ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericVector new_cent(ref_cent.length());
    Rcpp::IntegerVector num_vals(ref_cent.length());

    SEXP Y = Rcpp::List::create(Rcpp::_["cent"] = ref_cent);
    DtwBacktrackCalculator backtrack_calculator(DOTS, X, Y, is_multivariate);
    DbaUv parallel_worker(std::move(backtrack_calculator), new_cent, num_vals);
    int grain = get_grain(series.length(), num_threads);

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        // sum step
        RcppParallel::parallelFor(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = uv_average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
        Rcpp::checkUserInterrupt();
    }
    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'" << std::endl;
        Rflush();
    }
    return new_cent;
}

// =================================================================================================
/* multivariate DBA considering each variable separately */
// =================================================================================================

SEXP dba_mv_by_variable(const Rcpp::NumericMatrix& mv_ref_cent)
{
    ny = mv_ref_cent.nrow();
    nv = 1; // careful! this is used by the uv_* functions so it must be 1
    Rcpp::NumericMatrix mat_cent(ny, mv_ref_cent.ncol());
    Rcpp::NumericVector x(max_lengths_mv());
    Rcpp::NumericVector ref_cent(ny);
    Rcpp::NumericVector new_cent(ny);
    Rcpp::IntegerVector num_vals(ny);
    Rcpp::NumericVector kahan_c(ny), kahan_y(ny), kahan_t(ny); // for Kahan summation

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
            Rcpp::checkUserInterrupt();
        }

        if (iter > max_iter && trace) {
            Rcpp::Rcout << " Did not 'converge'" << std::endl;
            Rflush();
        }

        mat_cent(Rcpp::_, j) = new_cent;
    }

    return mat_cent;
}

// =================================================================================================
/* multivariate DBA considering each series as a whole */
// =================================================================================================

SEXP dba_mv_by_series(const Rcpp::NumericMatrix& centroid)
{
    Rcpp::NumericMatrix ref_cent = Rcpp::clone(centroid);
    ny = ref_cent.nrow();
    nv = ref_cent.ncol();
    Rcpp::NumericMatrix mat_cent(ny, nv);
    Rcpp::IntegerMatrix mat_vals(ny, nv);
    Rcpp::NumericMatrix kahan_c(ny, nv), kahan_y(ny, nv), kahan_t(ny, nv); // for Kahan summation

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        reset_matrices(mat_cent, mat_vals, kahan_c);

        // sum step
        for (int i = 0; i < series.length(); i++) {
            Rcpp::NumericMatrix x = series[i];
            nx = x.nrow();
            mv_set_alignment(x, ref_cent);
            mv_sum_step(x, mat_cent, mat_vals, kahan_c, kahan_y, kahan_t);
        }

        // average step with check for 'convergence' and update ref_cent
        bool converged = mv_average_step(mat_cent, mat_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
        Rcpp::checkUserInterrupt();
    }

    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'" << std::endl;
        Rflush();
    }

    return mat_cent;
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP dba(SEXP X, SEXP CENT,
                    SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
                    SEXP MV, SEXP MV_VER, SEXP DOTS, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    bool is_multivariate = Rcpp::as<bool>(MV);
    series = Rcpp::List(X);

    max_iter = Rcpp::as<int>(MAX_ITER);
    delta = Rcpp::as<double>(DELTA);
    trace = Rcpp::as<bool>(TRACE);
    num_threads = Rcpp::as<bool>(NUM_THREADS);

    Rcpp::List dots(DOTS);
    window = dots["window.size"];
    norm = dots["norm"];
    step = dots["step.pattern"];
    backtrack = dots["backtrack"];
    normalize = dots["normalize"];
    gcm = dots["gcm"];

    if (is_multivariate) {
        if (Rcpp::as<int>(MV_VER) == 1)
            return dba_mv_by_variable(CENT);
        else
            return dba_mv_by_series(CENT);
    }
    else {
        return dba_uv(X, CENT, DOTS, is_multivariate);
    }
    END_RCPP
}

} // namespace dtwclust
