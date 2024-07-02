#include "R-gateways.h"

#include <algorithm> // std::fill
#include <math.h> // exp
#include <utility> // std::move

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/calculators.h"
#include "../distances/details.h" // sdtw
#include "../utils/KahanSummer.h"
#include "../utils/ParallelWorker.h"
#include "../utils/SurrogateMatrix.h"
#include "../utils/TSTSList.h"
#include "../utils/utils.h" // get_grain, id_t

namespace dtwclust {

// =================================================================================================
/* helpers */
// =================================================================================================

void init_matrices(const id_t m, const id_t n,
                   SurrogateMatrix<double>& cm,
                   SurrogateMatrix<double>& dm,
                   SurrogateMatrix<double>& em)
{
    for (id_t i = 1; i <= m; i++) {
        dm(i-1, n) = 0;
        cm(i, n+1) = R_NegInf;
    }
    for (id_t j = 1; j <= n; j++) {
        dm(m, j-1) = 0;
        cm(m+1, j) = R_NegInf;
    }
    cm(m+1, n+1) = cm(m,n);
    dm(m,n) = 0;
    em.fill(0);
    em((m+1)%2, n+1) = 1;
}

void update_em(const id_t i, const id_t n, const double gamma,
               SurrogateMatrix<double>& cm,
               SurrogateMatrix<double>& dm,
               SurrogateMatrix<double>& em)
{
    double a, b, c;
    for (id_t j = n; j > 0; j--) {
        a = exp((cm(i+1, j) - cm(i,j) - dm(i, j-1)) / gamma);
        b = exp((cm(i, j+1) - cm(i,j) - dm(i-1, j)) / gamma);
        c = exp((cm(i+1, j+1) - cm(i,j) - dm(i,j)) / gamma);
        em(i%2, j) = a * em((i+1)%2, j) + b * em(i%2, j+1) + c * em((i+1)%2, j+1);
    }
}

// =================================================================================================
/* thread-safe soft-DTW calculator */
// =================================================================================================

class SdtwCentCalculator : public DistanceCalculator
{
public:
    // constructor
    SdtwCentCalculator(const Rcpp::List& x, const Rcpp::List& y, const double gamma)
        : DistanceCalculator("SDTW_CENT")
        , gamma_(gamma)
        , x_(x)
        , y_(y)
    {
        // set values of max_len_*_
        max_len_x_ = this->maxLength(x_);
        max_len_y_ = this->maxLength(y_);
    }

    // calculate for given indices
    double calculate(const id_t i, const id_t j) override {
        return this->calculate(x_[i], y_[j]);
    }

    // clone that sets helper matrices
    SdtwCentCalculator* clone() const override {
        SdtwCentCalculator* ptr = new SdtwCentCalculator(*this);
        ptr->cm_ = SurrogateMatrix<double>(max_len_x_ + 2, max_len_y_ + 2);
        ptr->dm_ = SurrogateMatrix<double>(max_len_x_ + 1, max_len_y_ + 1);
        return ptr;
    }

    // calculate
    double calculate(const arma::mat& x, const arma::mat& y) {
        if (!cm_ || !dm_) return -1;

        SurrogateMatrix<const double> temp_x(x.n_rows, x.n_cols, &x[0]);
        SurrogateMatrix<const double> temp_y(y.n_rows, y.n_cols, &y[0]);
        return sdtw(temp_x, temp_y, gamma_, cm_, dm_);
    }

    // input parameters
    double gamma_;
    // input series
    TSTSList<arma::mat> x_, y_;
    // helper "matrices"
    SurrogateMatrix<double> cm_, dm_;
    // to dimension cm_ and dm_
    int max_len_x_, max_len_y_;
};

// =================================================================================================
/* the parallel worker for the univariate version */
// =================================================================================================

class SdtwUv : public ParallelWorker {
public:
    // constructor
    SdtwUv(const SdtwCentCalculator&& dist_calculator,
           const Rcpp::NumericVector& weights,
           Rcpp::NumericVector& gradient,
           double& objective,
           const double gamma,
           const int grain)
        : ParallelWorker(grain, 1, 10)
        , dist_calculator_(dist_calculator)
        , weights_(weights)
        , gamma_(gamma)
        , gradient_summer_(&gradient[0], gradient.length())
        , objective_summer_(&objective, 1)
    { }

    // parallel loop across specified range
    void work_it(std::size_t begin, std::size_t end) override {
        // local copy of calculator so it is setup separately for each thread
        mutex_.lock();
        SdtwCentCalculator* local_calculator = dist_calculator_.clone();
        SurrogateMatrix<double> em(2, local_calculator->max_len_y_ + 2);
        mutex_.unlock();

        // calculations for these series
        //   Includes the calculation of soft-DTW gradient + jacobian product.
        SurrogateMatrix<double>& cm = local_calculator->cm_;
        SurrogateMatrix<double>& dm = local_calculator->dm_;
        const auto& x = local_calculator->x_[0];
        id_t m = x.n_rows;
        for (std::size_t id = begin; id < end; id++) {
            if (is_interrupted(id)) break; // nocov

            const auto& y = local_calculator->y_[id];
            double dist = local_calculator->calculate(0,id);

            mutex_.lock();
            objective_summer_.add(weights_[id] * dist, 0);
            mutex_.unlock();

            id_t n = y.n_rows;
            init_matrices(m, n, cm, dm, em);
            for (id_t i = m; i > 0; i--) {
                update_em(i, n, gamma_, cm, dm, em);
                double grad = 0;
                for (id_t j = 0; j < n; j++) grad += em(i%2, j+1) * 2 * (x[i-1] - y[j]);

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
};

// =================================================================================================
/* the parallel worker for the multivariate version */
// =================================================================================================

class SdtwMv : public ParallelWorker {
public:
    // constructor
    SdtwMv(const SdtwCentCalculator&& dist_calculator,
           const Rcpp::NumericVector& weights,
           Rcpp::NumericMatrix& gradient,
           double& objective,
           const double gamma,
           const int grain)
        : ParallelWorker(grain, 1, 10)
        , dist_calculator_(dist_calculator)
        , weights_(weights)
        , gamma_(gamma)
        , gradient_summer_(&gradient[0], gradient.nrow(), gradient.ncol())
        , objective_summer_(&objective, 1)
    { }

    // parallel loop across specified range
    void work_it(std::size_t begin, std::size_t end) override {
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
        const auto& x = local_calculator->x_[0];
        id_t m = x.n_rows, dim = x.n_cols;
        for (std::size_t id = begin; id < end; id++) {
            if (is_interrupted(id)) break; // nocov

            const auto& y = local_calculator->y_[id];
            double dist = local_calculator->calculate(0,id);

            mutex_.lock();
            objective_summer_.add(weights_[id] * dist, 0);
            if (!grad) grad = new double[dim];
            mutex_.unlock();

            id_t n = y.n_rows;
            init_matrices(m, n, cm, dm, em);
            for (id_t i = m; i > 0; i--) {
                update_em(i, n, gamma_, cm, dm, em);
                std::fill(grad, grad + dim, 0);
                for (id_t j = 0; j < n; j++) {
                    for (id_t k = 0; k < dim; k++) {
                        grad[k] += em(i%2, j+1) * 2 * (x.at(i-1, k) - y.at(j,k));
                    }
                }

                mutex_.lock();
                for (id_t k = 0; k < dim; k++)
                    gradient_summer_.add(weights_[id] * grad[k], i-1, k);
                mutex_.unlock();

                if (i == m) em((m+1)%2, n+1) = 0;
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
    int grain = get_grain(series.length(), num_threads);
    if (grain == DTWCLUST_MIN_GRAIN) grain = 1;

    // compute objective and gradient
    if (Rcpp::as<bool>(MV)) {
        Rcpp::NumericMatrix cent(CENTROID);
        Rcpp::NumericMatrix gradient(cent.nrow(), cent.ncol());

        double objective = 0;
        Rcpp::List x = Rcpp::List::create(CENTROID);

        SdtwCentCalculator dist_calculator(x, series, gamma);
        SdtwMv parallel_worker(std::move(dist_calculator), WEIGHTS, gradient, objective, gamma, grain);
        parallel_for(0, series.length(), parallel_worker, grain);

        return Rcpp::List::create(
            Rcpp::_["objective"] = objective,
            Rcpp::_["gradient"] = gradient
        );
    }
    else {
        Rcpp::NumericVector cent(CENTROID);
        Rcpp::NumericVector gradient(cent.length());

        double objective = 0;
        Rcpp::List x = Rcpp::List::create(CENTROID);

        SdtwCentCalculator dist_calculator(x, series, gamma);
        SdtwUv parallel_worker(std::move(dist_calculator), WEIGHTS, gradient, objective, gamma, grain);
        parallel_for(0, series.length(), parallel_worker, grain);

        return Rcpp::List::create(
            Rcpp::_["objective"] = objective,
            Rcpp::_["gradient"] = gradient
        );
    }
    END_RCPP
}

} // namespace dtwclust
