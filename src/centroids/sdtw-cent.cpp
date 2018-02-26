#include "centroids.h"

#include <algorithm> // std::fill
#include <math.h> // exp
#include <utility> // std::move

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distances/distances-details.h" // sdtw
#include "../utils/SurrogateMatrix.h"
#include "../utils/TSTSList.h"
#include "../utils/utils.h" // KahanSummer

namespace dtwclust {

// =================================================================================================
/* helpers */
// =================================================================================================

void init_matrices(const int m, const int n,
                   SurrogateMatrix<double>& cm,
                   SurrogateMatrix<double>& dm,
                   SurrogateMatrix<double>& em)
{
    for (int i = 1; i <= m; i++) {
        dm(i-1,n) = 0;
        cm(i,n+1) = R_NegInf;
    }
    for (int j = 1; j <= n; j++) {
        dm(m,j-1) = 0;
        cm(m+1,j) = R_NegInf;
    }
    cm(m+1,n+1) = cm(m,n);
    dm(m,n) = 0;
    em.fill(0);
    em((m+1)%2,n+1) = 1;
}

void update_em(const int i, const int n, const double gamma,
               SurrogateMatrix<double>& cm,
               SurrogateMatrix<double>& dm,
               SurrogateMatrix<double>& em)
{
    double a, b, c;
    for (int j = n; j > 0; j--) {
        a = exp((cm(i+1,j) - cm(i,j) - dm(i,j-1)) / gamma);
        b = exp((cm(i,j+1) - cm(i,j) - dm(i-1,j)) / gamma);
        c = exp((cm(i+1,j+1) - cm(i,j) - dm(i,j)) / gamma);
        em(i%2,j) = a * em((i+1)%2,j) + b * em(i%2,j+1) + c * em((i+1)%2,j+1);
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
        ptr->cm_ = SurrogateMatrix<double>(max_len_x_ + 2, max_len_y_ + 2);
        ptr->dm_ = SurrogateMatrix<double>(max_len_x_ + 1, max_len_y_ + 1);
        return ptr;
    }

    // univariate calculate
    double calculate(const RcppParallel::RVector<double>& x, const RcppParallel::RVector<double>& y) {
        if (!cm_ || !dm_) return -1;
        int nx = x.length();
        int ny = y.length();
        return sdtw(&x[0], &y[0], nx, ny, 1, gamma_, cm_, dm_);
    }

    // multivariate calculate
    double calculate(const RcppParallel::RMatrix<double>& x, const RcppParallel::RMatrix<double>& y) {
        if (!cm_ || !dm_) return -1;
        int nx = x.nrow();
        int ny = y.nrow();
        int num_var = x.ncol();
        return sdtw(&x[0], &y[0], nx, ny, num_var, gamma_, cm_, dm_);
    }

    // input parameters
    double gamma_;
    bool is_multivariate_;
    // input series (univariate)
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    // input series (multivariate)
    TSTSList<Rcpp::NumericMatrix> x_mv_, y_mv_;
    // helper "matrices"
    SurrogateMatrix<double> cm_, dm_;
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
        SurrogateMatrix<double> em(2, local_calculator->max_len_y_ + 2);
        mutex_.unlock();
        // calculations for these series
        //   Includes the calculation of soft-DTW gradient + jacobian product.
        SurrogateMatrix<double>& cm = local_calculator->cm_;
        SurrogateMatrix<double>& dm = local_calculator->dm_;
        const RcppParallel::RVector<double>& x = local_calculator->x_uv_[0];
        int m = x.length();
        for (std::size_t id = begin; id < end; id++) {
            const RcppParallel::RVector<double>& y = local_calculator->y_uv_[id];
            double dist = local_calculator->calculate(0,id);

            mutex_.lock();
            objective_summer_.add(weights_[id] * dist, 0);
            mutex_.unlock();

            int n = y.length();
            init_matrices(m, n, cm, dm, em);
            for (int i = m; i > 0; i--) {
                update_em(i, n, gamma_, cm, dm, em);
                double grad = 0;
                for (int j = 0; j < n; j++) grad += em(i%2,j+1) * 2 * (x[i-1] - y[j]);

                mutex_.lock();
                gradient_summer_.add(weights_[id] * grad, i-1);
                mutex_.unlock();

                if (i == m) em((m+1)%2,n+1) = 0;
            }
        }
        // finish
        mutex_.lock();
        delete local_calculator;
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
        SurrogateMatrix<double> em(2, local_calculator->max_len_y_ + 2);
        mutex_.unlock();
        // calculations for these series
        //   Includes the calculation of soft-DTW gradient + jacobian product.
        double* grad = nullptr;
        SurrogateMatrix<double>& cm = local_calculator->cm_;
        SurrogateMatrix<double>& dm = local_calculator->dm_;
        const RcppParallel::RMatrix<double>& x = local_calculator->x_mv_[0];
        int m = x.nrow(), dim = x.ncol();
        for (std::size_t id = begin; id < end; id++) {
            const RcppParallel::RMatrix<double>& y = local_calculator->y_mv_[id];
            double dist = local_calculator->calculate(0,id);

            mutex_.lock();
            objective_summer_.add(weights_[id] * dist, 0);
            if (!grad) grad = new double[dim];
            mutex_.unlock();

            int n = y.nrow();
            init_matrices(m, n, cm, dm, em);
            for (int i = m; i > 0; i--) {
                update_em(i, n, gamma_, cm, dm, em);
                std::fill(grad, grad + dim, 0);
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < dim; k++) {
                        grad[k] += em(i%2,j+1) * 2 * (x(i-1,k) - y(j,k));
                    }
                }

                mutex_.lock();
                for (int k = 0; k < dim; k++) gradient_summer_.add(weights_[id] * grad[k], i-1, k);
                mutex_.unlock();

                if (i == m) em((m+1)%2,n+1) = 0;
            }
        }
        // finish
        mutex_.lock();
        delete local_calculator;
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

extern "C" SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID,
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
