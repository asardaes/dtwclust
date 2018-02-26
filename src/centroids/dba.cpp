#include "centroids.h"

#include <cmath> // std::abs
#include <cstddef> // std::size_t
#include <utility> // std::move

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distances/distances-details.h" // dtw_basic_par
#include "../utils/TSTSList.h"
#include "../utils/utils.h" // KahanSummer, Rflush

namespace dtwclust {

// =================================================================================================
/* shared variables */
// =================================================================================================

static int max_iter, num_threads;
static double delta;
static bool trace;

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

    // default copy/move constructor (must be explicitly defined due to custom destructor)
    DtwBacktrackCalculator(const DtwBacktrackCalculator&) = default;
    DtwBacktrackCalculator(DtwBacktrackCalculator&&) = default;

    // calculate for given indices (inherited)
    double calculate(const int i, const int j) override {
        if (is_multivariate_)
            return this->calculate(x_mv_[i], y_mv_[j]);
        else
            return this->calculate(x_uv_[i], y_uv_[j]);
    }

    // calculate for given indices (custom for multivariate by-variable version)
    double calculate(const int i, const int j, const int k) {
        return this->calculate(x_mv_[i], y_mv_[j], k);
    }

    // clone to setup helpers in each thread
    DtwBacktrackCalculator* clone() const override {
        DtwBacktrackCalculator* ptr = new DtwBacktrackCalculator(*this);
        ptr->gcm_ = new double[(max_len_x_ + 1) * (max_len_y_ + 1)];
        ptr->index1_ = new int[max_len_x_ + max_len_y_];
        ptr->index2_ = new int[max_len_x_ + max_len_y_];
        return ptr;
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
        return dtw_basic_par(&x[0], &y[0],
                             nx, ny, 1,
                             window_, norm_, step_, normalize_,
                             gcm_, true, index1_, index2_, &path_);
    }

    // by-series multivariate calculate
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

    // by-variable multivariate calculate
    double calculate(const RcppParallel::RMatrix<double>& x,
                     const RcppParallel::RMatrix<double>& y,
                     const int k)
    {
        if (!gcm_ || !index1_ || !index2_) return -1;
        int nx = x.nrow();
        int ny = y.nrow();
        return dtw_basic_par(&x[0] + (nx * k), &y[0] + (ny * k),
                             nx, ny, 1,
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
        , summer_(&new_cent_[0], new_cent_.length())
    { }

    // for reusability
    void reset() {
        summer_.reset();
    }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();
        // kahan sum step
        for (std::size_t i = begin; i < end; i++) {
            local_calculator->calculate(i,0);
            const RcppParallel::RVector<double>& x = local_calculator->x_uv_[i];
            mutex_.lock();
            for (int ii = local_calculator->path_ - 1; ii >= 0; ii--) {
                int i1 = local_calculator->index1_[ii] - 1;
                int i2 = local_calculator->index2_[ii] - 1;
                summer_.add(x[i1], i2);
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
    // sum helper
    KahanSummer summer_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// =================================================================================================
/* the parallel worker for the multivariate by-series version */
// =================================================================================================

class DbaMvBySeries : public RcppParallel::Worker {
public:
    // constructor
    DbaMvBySeries(const DtwBacktrackCalculator&& backtrack_calculator,
                  const Rcpp::NumericMatrix& new_cent,
                  const Rcpp::IntegerMatrix& num_vals)
        : backtrack_calculator_(backtrack_calculator)
        , new_cent_(new_cent)
        , num_vals_(num_vals)
        , summer_(&new_cent_[0], new_cent_.nrow(), new_cent_.ncol())
    { }

    // for reusability
    void reset() {
        summer_.reset();
    }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();
        // kahan sum step
        for (std::size_t i = begin; i < end; i++) {
            local_calculator->calculate(i,0);
            const RcppParallel::RMatrix<double>& x = local_calculator->x_mv_[i];
            mutex_.lock();
            for (int j = 0; j < static_cast<int>(new_cent_.ncol()); j++) {
                for (int ii = local_calculator->path_ - 1; ii >= 0; ii--) {
                    int i1 = local_calculator->index1_[ii] - 1;
                    int i2 = local_calculator->index2_[ii] - 1;
                    summer_.add(x(i1,j), i2, j);
                    num_vals_(i2,j) += 1;
                }
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
    RcppParallel::RMatrix<double> new_cent_;
    RcppParallel::RMatrix<int> num_vals_;
    // sum helper
    KahanSummer summer_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// =================================================================================================
/* the parallel worker for the multivariate by-variable version */
// =================================================================================================

class DbaMvByVariable : public RcppParallel::Worker {
public:
    // constructor
    DbaMvByVariable(const DtwBacktrackCalculator&& backtrack_calculator,
                    const Rcpp::NumericMatrix& new_cent,
                    const Rcpp::IntegerMatrix& num_vals)
        : backtrack_calculator_(backtrack_calculator)
        , new_cent_(new_cent)
        , num_vals_(num_vals)
        , summer_(&new_cent_[0], new_cent_.nrow(), new_cent_.ncol())
    { }

    // for reusability
    void reset() {
        summer_.reset();
    }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();
        // kahan sum step
        for (std::size_t i = begin; i < end; i++) {
            const RcppParallel::RMatrix<double>& x = local_calculator->x_mv_[i];
            for (int j = 0; j < static_cast<int>(new_cent_.ncol()); j++) {
                local_calculator->calculate(i,0,j);
                mutex_.lock();
                for (int ii = local_calculator->path_ - 1; ii >= 0; ii--) {
                    int i1 = local_calculator->index1_[ii] - 1;
                    int i2 = local_calculator->index2_[ii] - 1;
                    summer_.add(x(i1,j), i2, j);
                    num_vals_(i2,j) += 1;
                }
                mutex_.unlock();
            }
        }
        // finish
        mutex_.lock();
        delete local_calculator;
        mutex_.unlock();
    }

private:
    const DtwBacktrackCalculator& backtrack_calculator_;
    RcppParallel::RMatrix<double> new_cent_;
    RcppParallel::RMatrix<int> num_vals_;
    // sum helper
    KahanSummer summer_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// =================================================================================================
/* average step with check for 'convergence' and update ref_cent for vectors and matrices */
// =================================================================================================

// univariate
bool average_step(Rcpp::NumericVector& new_cent,
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

// multivariate
bool average_step(Rcpp::NumericMatrix& new_cent,
                  const Rcpp::IntegerMatrix& num_vals,
                  Rcpp::NumericMatrix& ref_cent)
{
    bool converged = true;
    for (int j = 0; j < new_cent.ncol(); j++) {
        for (int i = 0; i < new_cent.nrow(); i++) {
            new_cent(i,j) /= num_vals(i,j);
            if (std::abs(new_cent(i,j) - ref_cent(i,j)) >= delta) converged = false;
            ref_cent(i,j) = new_cent(i,j);
        }
    }
    return converged;
}

// =================================================================================================
/* tracing */
// =================================================================================================

void print_trace(const bool converged, const int iter)
{
    if (trace) {
        if (converged) {
            Rcpp::Rcout << " " << iter << " - Converged!" << std::endl;
        }
        else {
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

SEXP dba_uv(const SEXP& X, const Rcpp::NumericVector& centroid, const SEXP& DOTS)
{
    Rcpp::List series(X);
    Rcpp::NumericVector ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericVector new_cent(ref_cent.length());
    Rcpp::IntegerVector num_vals(ref_cent.length());

    SEXP Y = Rcpp::List::create(Rcpp::_["cent"] = ref_cent);
    DtwBacktrackCalculator backtrack_calculator(DOTS, X, Y, false);
    DbaUv parallel_worker(std::move(backtrack_calculator), new_cent, num_vals);
    int grain = series.length() / num_threads;
    if (grain == 0) grain = 1;

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        parallel_worker.reset();
        // sum step
        RcppParallel::parallelFor(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
        Rcpp::checkUserInterrupt();
    }
    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'\n";
        Rflush();
    }
    return new_cent;
}

// =================================================================================================
/* multivariate DBA considering each variable separately */
// =================================================================================================

SEXP dba_mv_by_variable(const SEXP& X, const Rcpp::NumericMatrix& centroid, const SEXP& DOTS)
{
    Rcpp::List series(X);
    Rcpp::NumericMatrix ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericMatrix new_cent(ref_cent.nrow(), ref_cent.ncol());
    Rcpp::IntegerMatrix num_vals(ref_cent.nrow(), ref_cent.ncol());

    SEXP Y = Rcpp::List::create(Rcpp::_["cent"] = ref_cent);
    DtwBacktrackCalculator backtrack_calculator(DOTS, X, Y, true);
    DbaMvByVariable parallel_worker(std::move(backtrack_calculator), new_cent, num_vals);
    int grain = series.length() / num_threads;
    if (grain == 0) grain = 1;

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        parallel_worker.reset();
        // sum step
        RcppParallel::parallelFor(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
        Rcpp::checkUserInterrupt();
    }
    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'\n";
        Rflush();
    }
    return new_cent;
}

// =================================================================================================
/* multivariate DBA considering each series as a whole */
// =================================================================================================

SEXP dba_mv_by_series(const SEXP& X, const Rcpp::NumericMatrix& centroid, const SEXP& DOTS)
{
    Rcpp::List series(X);
    Rcpp::NumericMatrix ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericMatrix new_cent(ref_cent.nrow(), ref_cent.ncol());
    Rcpp::IntegerMatrix num_vals(ref_cent.nrow(), ref_cent.ncol());

    SEXP Y = Rcpp::List::create(Rcpp::_["cent"] = ref_cent);
    DtwBacktrackCalculator backtrack_calculator(DOTS, X, Y, true);
    DbaMvBySeries parallel_worker(std::move(backtrack_calculator), new_cent, num_vals);
    int grain = series.length() / num_threads;
    if (grain == 0) grain = 1;

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        parallel_worker.reset();
        // sum step
        RcppParallel::parallelFor(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
        Rcpp::checkUserInterrupt();
    }
    if (iter > max_iter && trace) {
        Rcpp::Rcout << " Did not 'converge'\n";
        Rflush();
    }
    return new_cent;
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

extern "C" SEXP dba(SEXP X, SEXP CENT,
                    SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
                    SEXP MV, SEXP MV_VER, SEXP DOTS, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    max_iter = Rcpp::as<int>(MAX_ITER);
    delta = Rcpp::as<double>(DELTA);
    trace = Rcpp::as<bool>(TRACE);
    num_threads = Rcpp::as<int>(NUM_THREADS);
    if (Rcpp::as<bool>(MV)) {
        if (Rcpp::as<int>(MV_VER) == 1)
            return dba_mv_by_variable(X, CENT, DOTS);
        else
            return dba_mv_by_series(X, CENT, DOTS);
    }
    else {
        return dba_uv(X, CENT, DOTS);
    }
    END_RCPP
}

} // namespace dtwclust
