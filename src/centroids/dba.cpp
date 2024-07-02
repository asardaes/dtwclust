#include "R-gateways.h"

#include <cmath> // std::abs
#include <utility> // std::move

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/calculators.h"
#include "../distances/details.h" // dtw_basi
#include "../utils/KahanSummer.h"
#include "../utils/ParallelWorker.h"
#include "../utils/SurrogateMatrix.h"
#include "../utils/TSTSList.h"
#include "../utils/utils.h" // Rflush, get_grain, id_t

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
    DtwBacktrackCalculator(const Rcpp::List& dist_args, const Rcpp::List& x, const Rcpp::List& y)
        : DistanceCalculator("DTW_BACTRACK")
        , x_(x)
        , y_(y)
    {
        window_ = Rcpp::as<int>(dist_args["window.size"]);
        norm_ = Rcpp::as<double>(dist_args["norm"]);
        step_ = Rcpp::as<double>(dist_args["step.pattern"]);
        normalize_ = Rcpp::as<bool>(dist_args["normalize"]);
        // set value of max_len_*_
        max_len_x_ = this->maxLength(x_);
        max_len_y_ = this->maxLength(y_);
    }

    // calculate for given indices (inherited)
    double calculate(const id_t i, const id_t j) override {
        return this->calculate(x_[i], y_[j]);
    }

    // calculate for given indices (custom for multivariate by-variable version)
    double calculate(const id_t i, const id_t j, const id_t k) {
        return this->calculate(x_[i], y_[j], k);
    }

    // clone to setup helpers in each thread
    DtwBacktrackCalculator* clone() const override {
        DtwBacktrackCalculator* ptr = new DtwBacktrackCalculator(*this);
        ptr->lcm_ = SurrogateMatrix<double>(max_len_x_ + 1, max_len_y_ + 1);
        ptr->index1_ = SurrogateMatrix<int>(max_len_x_ + max_len_y_, 1);
        ptr->index2_ = SurrogateMatrix<int>(max_len_x_ + max_len_y_, 1);
        return ptr;
    }

    // input series
    TSTSList<arma::mat> x_, y_;
    // helpers for backtracking
    int path_;
    SurrogateMatrix<int> index1_, index2_;

private:
    // primary calculate
    double calculate(const arma::mat& x, const arma::mat& y) {
        if (!lcm_ || !index1_ || !index2_) return -1;

        SurrogateMatrix<const double> temp_x(x.n_rows, x.n_cols, &x[0]);
        SurrogateMatrix<const double> temp_y(y.n_rows, y.n_cols, &y[0]);
        return dtw_basic(lcm_, temp_x, temp_y,
                         window_, norm_, step_, normalize_, true,
                         index1_, index2_, path_);
    }

    // by-variable multivariate calculate
    double calculate(const arma::mat& x, const arma::mat& y, const id_t k) {
        if (!lcm_ || !index1_ || !index2_) return -1;

        SurrogateMatrix<const double> temp_x(x.n_rows, 1, &x[0] + (k * x.n_rows));
        SurrogateMatrix<const double> temp_y(y.n_rows, 1, &y[0] + (k * y.n_rows));
        return dtw_basic(lcm_, temp_x, temp_y,
                         window_, norm_, step_, normalize_, true,
                         index1_, index2_, path_);
    }

    // input parameters
    int window_;
    double norm_, step_;
    bool normalize_;
    // helper "matrix"
    SurrogateMatrix<double> lcm_;
    // to dimension helpers
    int max_len_x_, max_len_y_;
};

// =================================================================================================
/* the parallel worker for the univariate version */
// =================================================================================================

class DbaUv : public ParallelWorker {
public:
    // constructor
    DbaUv(const DtwBacktrackCalculator&& backtrack_calculator,
          const Rcpp::NumericVector& new_cent,
          const Rcpp::IntegerVector& num_vals,
          const int grain)
        : ParallelWorker(grain, 50, 100)
        , backtrack_calculator_(backtrack_calculator)
        , new_cent_(new_cent)
        , num_vals_(num_vals)
        , summer_(&new_cent_[0], new_cent_.length())
    { }

    // for reusability
    void reset() {
        summer_.reset();
    }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();

        // kahan sum step
        for (id_t i = begin; i < end; i++) {
            if (is_interrupted(i)) break; // nocov

            local_calculator->calculate(i,0);
            const auto& x = local_calculator->x_[i];

            mutex_.lock();
            for (int ii = 0; ii < local_calculator->path_; ii++) {
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
};

// =================================================================================================
/* the parallel worker for the multivariate by-series version */
// =================================================================================================

class DbaMvBySeries : public ParallelWorker {
public:
    // constructor
    DbaMvBySeries(const DtwBacktrackCalculator&& backtrack_calculator,
                  const Rcpp::NumericMatrix& new_cent,
                  const Rcpp::IntegerMatrix& num_vals,
                  const int grain)
        : ParallelWorker(grain, 50, 100)
        , backtrack_calculator_(backtrack_calculator)
        , new_cent_(new_cent)
        , num_vals_(num_vals)
        , summer_(&new_cent_[0], new_cent_.nrow(), new_cent_.ncol())
    { }

    // for reusability
    void reset() {
        summer_.reset();
    }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();

        // kahan sum step
        for (id_t i = begin; i < end; i++) {
            if (is_interrupted(i)) break; // nocov

            local_calculator->calculate(i,0);
            const auto& x = local_calculator->x_[i];

            mutex_.lock();
            for (id_t j = 0; j < new_cent_.ncol(); j++) {
                for (int ii = 0; ii < local_calculator->path_; ii++) {
                    int i1 = local_calculator->index1_[ii] - 1;
                    int i2 = local_calculator->index2_[ii] - 1;
                    summer_.add(x.at(i1,j), i2, j);
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
};

// =================================================================================================
/* the parallel worker for the multivariate by-variable version */
// =================================================================================================

class DbaMvByVariable : public ParallelWorker {
public:
    // constructor
    DbaMvByVariable(const DtwBacktrackCalculator&& backtrack_calculator,
                    const Rcpp::NumericMatrix& new_cent,
                    const Rcpp::IntegerMatrix& num_vals,
                    const int grain)
        : ParallelWorker(grain, 10, 50)
        , backtrack_calculator_(backtrack_calculator)
        , new_cent_(new_cent)
        , num_vals_(num_vals)
        , summer_(&new_cent_[0], new_cent_.nrow(), new_cent_.ncol())
    { }

    // for reusability
    void reset() {
        summer_.reset();
    }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        DtwBacktrackCalculator* local_calculator = backtrack_calculator_.clone();
        mutex_.unlock();

        // kahan sum step
        for (id_t i = begin; i < end; i++) {
            if (is_interrupted(i)) break; // nocov

            const auto& x = local_calculator->x_[i];
            for (id_t j = 0; j < new_cent_.ncol(); j++) {
                local_calculator->calculate(i,0,j);

                mutex_.lock();
                for (int ii = 0; ii < local_calculator->path_; ii++) {
                    int i1 = local_calculator->index1_[ii] - 1;
                    int i2 = local_calculator->index2_[ii] - 1;
                    summer_.add(x.at(i1,j), i2, j);
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
};

// =================================================================================================
/* average step with check for 'convergence' and update ref_cent for vectors and matrices */
// =================================================================================================

// univariate
bool average_step(Rcpp::NumericVector& new_cent,
                  const Rcpp::IntegerVector& num_vals,
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

SEXP dba_uv(const Rcpp::List& series, const Rcpp::NumericVector& centroid, const SEXP& DOTS)
{
    Rcpp::NumericVector ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericVector new_cent(ref_cent.length());
    Rcpp::IntegerVector num_vals(ref_cent.length());

    DtwBacktrackCalculator backtrack_calculator(DOTS, series, Rcpp::List::create(ref_cent));
    int grain = get_grain(series.length(), num_threads);
    if (grain == DTWCLUST_MIN_GRAIN) grain = 1;
    DbaUv parallel_worker(std::move(backtrack_calculator), new_cent, num_vals, grain);

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        parallel_worker.reset();
        // sum step
        parallel_for(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
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

SEXP dba_mv_by_variable(const Rcpp::List& series, const Rcpp::NumericMatrix& centroid, const SEXP& DOTS)
{
    Rcpp::NumericMatrix ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericMatrix new_cent(ref_cent.nrow(), ref_cent.ncol());
    Rcpp::IntegerMatrix num_vals(ref_cent.nrow(), ref_cent.ncol());

    DtwBacktrackCalculator backtrack_calculator(DOTS, series, Rcpp::List::create(ref_cent));
    int grain = get_grain(series.length(), num_threads);
    if (grain == DTWCLUST_MIN_GRAIN) grain = 1;
    DbaMvByVariable parallel_worker(std::move(backtrack_calculator), new_cent, num_vals, grain);

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        parallel_worker.reset();
        // sum step
        parallel_for(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
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

SEXP dba_mv_by_series(const Rcpp::List& series, const Rcpp::NumericMatrix& centroid, const SEXP& DOTS)
{
    Rcpp::NumericMatrix ref_cent = Rcpp::clone(centroid);
    Rcpp::NumericMatrix new_cent(ref_cent.nrow(), ref_cent.ncol());
    Rcpp::IntegerMatrix num_vals(ref_cent.nrow(), ref_cent.ncol());

    DtwBacktrackCalculator backtrack_calculator(DOTS, series, Rcpp::List::create(ref_cent));
    int grain = get_grain(series.length(), num_threads);
    if (grain == DTWCLUST_MIN_GRAIN) grain = 1;
    DbaMvBySeries parallel_worker(std::move(backtrack_calculator), new_cent, num_vals, grain);

    if (trace) Rcpp::Rcout << "\tDBA Iteration:";
    int iter = 1;
    while (iter <= max_iter) {
        new_cent.fill(0);
        num_vals.fill(0);
        parallel_worker.reset();
        // sum step
        parallel_for(0, series.length(), parallel_worker, grain);
        // average step with check for 'convergence' and update ref_cent
        bool converged = average_step(new_cent, num_vals, ref_cent);
        print_trace(converged, iter);
        if (converged) break;
        iter++;
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
