#include "centroids.h"

#include <algorithm> // std::fill
#include <math.h> // exp

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distances/distances.h" // sdtw
#include "../utils/TSTSList.h"
#include "../utils/utils.h" // KahanSummer, d2s

namespace dtwclust {

// =================================================================================================
/* helpers */
// =================================================================================================

void init_matrices(const int m, const int n,
                   const int cm_nrows, const int dm_nrows, const int max_len_y,
                   double * const cm, double * const dm, double * const em)
{
    for (int i = 1; i <= m; i++) {
        dm[d2s(i-1, n, dm_nrows)] = 0;
        cm[d2s(i, n+1, cm_nrows)] = R_NegInf;
    }
    for (int j = 1; j <= n; j++) {
        dm[d2s(m, j-1, dm_nrows)] = 0;
        cm[d2s(m+1, j, cm_nrows)] = R_NegInf;
    }
    cm[d2s(m+1, n+1, cm_nrows)] = cm[d2s(m,n,cm_nrows)];
    dm[d2s(m,n,dm_nrows)] = 0;
    std::fill(em, em + 2 * (max_len_y + 2), 0);
    em[d2s((m+1)%2, n+1, 2)] = 1;
}

void update_em(const int i, const int n, const double gamma,
               const int cm_nrows, const int dm_nrows,
               const double * const cm, const double * const dm, double * const em)
{
    double a, b, c;
    for (int j = n; j > 0; j--) {
        a = exp((cm[d2s(i+1,j,cm_nrows)] - cm[d2s(i,j,cm_nrows)] - dm[d2s(i,j-1,dm_nrows)]) / gamma);
        b = exp((cm[d2s(i,j+1,cm_nrows)] - cm[d2s(i,j,cm_nrows)] - dm[d2s(i-1,j,dm_nrows)]) / gamma);
        c = exp((cm[d2s(i+1,j+1,cm_nrows)] - cm[d2s(i,j,cm_nrows)] - dm[d2s(i,j,dm_nrows)]) / gamma);
        em[d2s(i%2,j,2)] = a * em[d2s((i+1)%2,j,2)] + b * em[d2s(i%2,j+1,2)] + c * em[d2s((i+1)%2,j+1,2)];
    }
}

// =================================================================================================
/* thread-safe soft-DTW calculator */
// =================================================================================================

class SdtwCentCalculator : public DistanceCalculator
{
public:
    // constructor
    SdtwCentCalculator(const Rcpp::List& x, const Rcpp::List& y, const double gamma, const bool mv)
        : gamma_(gamma)
        , is_multivariate_(mv)
    {
        if (is_multivariate_) {
            x_mv_ = TSTSList<Rcpp::NumericMatrix>(x);
            y_mv_ = TSTSList<Rcpp::NumericMatrix>(y);
        }
        else {
            x_uv_ = TSTSList<Rcpp::NumericVector>(x);
            y_uv_ = TSTSList<Rcpp::NumericVector>(y);
        }
        // set values of max_len_*_
        max_len_x_ = this->maxLength(x, is_multivariate_);
        max_len_y_ = this->maxLength(y, is_multivariate_);
        // make sure pointers are null
        cm_ = nullptr;
        dm_ = nullptr;
    }

    // destructor
    ~SdtwCentCalculator() {
        if (cm_) delete[] cm_;
        if (dm_) delete[] dm_;
    }

    // calculate for given indices
    double calculate(const int i, const int j) override {
        if (is_multivariate_)
            return this->calculate(x_mv_[i], y_mv_[j]);
        else
            return this->calculate(x_uv_[i], y_uv_[j]);
    }

    // clone that sets helper matrices
    SdtwCentCalculator* clone() const override {
        SdtwCentCalculator* ptr = new SdtwCentCalculator(*this);
        ptr->cm_ = new double[(max_len_x_ + 2) * (max_len_y_ + 2)];
        ptr->dm_ = new double[(max_len_x_ + 1) * (max_len_y_ + 1)];
        return ptr;
    }

    // limits
    int xLimit() const override { return is_multivariate_ ? x_mv_.length() : x_uv_.length(); } // nocov start
    int yLimit() const override { return is_multivariate_ ? y_mv_.length() : y_uv_.length(); } // nocov end

    // univariate calculate
    double calculate(const RcppParallel::RVector<double>& x, const RcppParallel::RVector<double>& y) {
        if (!cm_ || !dm_) return -1;
        int nx = x.length();
        int ny = y.length();
        return sdtw(&x[0], &y[0], nx, ny, 1, gamma_, cm_, true, dm_);
    }

    // multivariate calculate
    double calculate(const RcppParallel::RMatrix<double>& x, const RcppParallel::RMatrix<double>& y) {
        if (!cm_ || !dm_) return -1;
        int nx = x.nrow();
        int ny = y.nrow();
        int num_var = x.ncol();
        return sdtw(&x[0], &y[0], nx, ny, num_var, gamma_, cm_, true, dm_);
    }

    // input parameters
    double gamma_;
    bool is_multivariate_;
    // input series (univariate)
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    // input series (multivariate)
    TSTSList<Rcpp::NumericMatrix> x_mv_, y_mv_;
    // helper "matrices"
    double *cm_, *dm_;
    // to dimension cm_ and dm_
    int max_len_x_, max_len_y_;
};

// =================================================================================================
/* the parallel worker for the univariate version */
// =================================================================================================

class SdtwUv : public RcppParallel::Worker {
public:
    // constructor
    SdtwUv(const SdtwCentCalculator&& dist_calculator,
           const Rcpp::NumericVector& weights,
           Rcpp::NumericVector& gradient,
           double& objective,
           const double gamma)
        : dist_calculator_(dist_calculator)
        , weights_(weights)
        , gamma_(gamma)
        , gradient_summer_(&gradient[0], gradient.length())
        , objective_summer_(&objective, 1)
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        SdtwCentCalculator* local_calculator = dist_calculator_.clone();
        double* em = new double[2 * (local_calculator->max_len_y_ + 2)];
        mutex_.unlock();
        // calculations for these series
        //   Includes the calculation of soft-DTW gradient + jacobian product.
        double* cm = local_calculator->cm_;
        double* dm = local_calculator->dm_;
        int cm_nrows = local_calculator->max_len_x_ + 2;
        int dm_nrows = local_calculator->max_len_x_ + 1;
        const RcppParallel::RVector<double>& x = local_calculator->x_uv_[0];
        int m = x.length();
        for (std::size_t id = begin; id < end; id++) {
            const RcppParallel::RVector<double>& y = local_calculator->y_uv_[id];
            double dist = local_calculator->calculate(0,id);
            mutex_.lock();
            objective_summer_.add(weights_[id] * dist, 0);
            mutex_.unlock();
            int n = y.length();
            init_matrices(m, n, cm_nrows, dm_nrows, local_calculator->max_len_y_, cm, dm, em);
            for (int i = m; i > 0; i--) {
                update_em(i, n, gamma_, cm_nrows, dm_nrows, cm, dm, em);
                double grad = 0;
                for (int j = 0; j < n; j++) grad += em[d2s(i%2, j+1, 2)] * 2 * (x[i-1] - y[j]);
                mutex_.lock();
                gradient_summer_.add(weights_[id] * grad, i-1);
                mutex_.unlock();
                if (i == m) em[d2s((m+1)%2, n+1, 2)] = 0;
            }
        }
        // finish
        mutex_.lock();
        delete local_calculator;
        delete[] em;
        mutex_.unlock();
    }

private:
    const SdtwCentCalculator& dist_calculator_;
    const RcppParallel::RVector<double> weights_;
    double gamma_;
    KahanSummer gradient_summer_, objective_summer_;
    tthread::mutex mutex_;
};

// =================================================================================================
/* the parallel worker for the multivariate version */
// =================================================================================================

class SdtwMv : public RcppParallel::Worker {
public:
    // constructor
    SdtwMv(const SdtwCentCalculator&& dist_calculator,
           const Rcpp::NumericVector& weights,
           Rcpp::NumericMatrix& gradient,
           double& objective,
           const double gamma)
        : dist_calculator_(dist_calculator)
        , weights_(weights)
        , gamma_(gamma)
        , gradient_summer_(&gradient[0], gradient.nrow(), gradient.ncol())
        , objective_summer_(&objective, 1)
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        SdtwCentCalculator* local_calculator = dist_calculator_.clone();
        double* em = new double[2 * (local_calculator->max_len_y_ + 2)];
        mutex_.unlock();
        // calculations for these series
        //   Includes the calculation of soft-DTW gradient + jacobian product.
        double* grad = nullptr;
        double* cm = local_calculator->cm_;
        double* dm = local_calculator->dm_;
        int cm_nrows = local_calculator->max_len_x_ + 2;
        int dm_nrows = local_calculator->max_len_x_ + 1;
        const RcppParallel::RMatrix<double>& x = local_calculator->x_mv_[0];
        int m = x.nrow(), dim = x.ncol();
        for (std::size_t id = begin; id < end; id++) {
            const RcppParallel::RMatrix<double>& y = local_calculator->y_mv_[id];
            double dist = local_calculator->calculate(0,id);
            mutex_.lock();
            objective_summer_.add(weights_[id] * dist, 0);
            mutex_.unlock();
            if (!grad) {
                mutex_.lock();
                grad = new double[dim];
                mutex_.unlock();
            }
            int n = y.nrow();
            init_matrices(m, n, cm_nrows, dm_nrows, local_calculator->max_len_y_, cm, dm, em);
            for (int i = m; i > 0; i--) {
                update_em(i, n, gamma_, cm_nrows, dm_nrows, cm, dm, em);
                std::fill(grad, grad + dim, 0);
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < dim; k++) {
                        grad[k] += em[d2s(i%2, j+1, 2)] * 2 * (x[d2s(i-1,k,m)] - y[d2s(j,k,n)]);
                    }
                }
                mutex_.lock();
                for (int k = 0; k < dim; k++) gradient_summer_.add(weights_[id] * grad[k], i-1, k);
                mutex_.unlock();
                if (i == m) em[d2s((m+1)%2, n+1, 2)] = 0;
            }
        }
        // finish
        mutex_.lock();
        delete local_calculator;
        delete[] em;
        if (grad) delete[] grad;
        mutex_.unlock();
    }

private:
    const SdtwCentCalculator& dist_calculator_;
    const RcppParallel::RVector<double> weights_;
    double gamma_;
    KahanSummer gradient_summer_, objective_summer_;
    tthread::mutex mutex_;
};

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID,
                          SEXP GAMMA, SEXP WEIGHTS,
                          SEXP MV, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    Rcpp::List series(SERIES);
    double gamma = Rcpp::as<double>(GAMMA);
    int num_threads = Rcpp::as<int>(NUM_THREADS);
    int grain = series.length() / num_threads;
    if (grain == 0) grain = 1;
    // compute objective and gradient
    if (Rcpp::as<bool>(MV)) {
        Rcpp::NumericMatrix cent(CENTROID);
        Rcpp::NumericMatrix gradient(cent.nrow(), cent.ncol());
        double objective = 0;
        Rcpp::List x = Rcpp::List::create(Rcpp::_["cent"] = CENTROID);
        SdtwCentCalculator dist_calculator(x, series, gamma, true);
        SdtwMv parallel_worker(std::move(dist_calculator), WEIGHTS, gradient, objective, gamma);
        RcppParallel::parallelFor(0, series.length(), parallel_worker, grain);
        return Rcpp::List::create(
            Rcpp::_["objective"] = objective,
            Rcpp::_["gradient"] = gradient
        );
    }
    else {
        Rcpp::NumericVector cent(CENTROID);
        Rcpp::NumericVector gradient(cent.length());
        double objective = 0;
        Rcpp::List x = Rcpp::List::create(Rcpp::_["cent"] = CENTROID);
        SdtwCentCalculator dist_calculator(x, series, gamma, false);
        SdtwUv parallel_worker(std::move(dist_calculator), WEIGHTS, gradient, objective, gamma);
        RcppParallel::parallelFor(0, series.length(), parallel_worker, grain);
        return Rcpp::List::create(
            Rcpp::_["objective"] = objective,
            Rcpp::_["gradient"] = gradient
        );
    }
    END_RCPP
}

} // namespace dtwclust
