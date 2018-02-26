#include "distmat-loops.h"

#include <cstddef> // std::size_t
#include <memory> // shared_ptr

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distance-calculators/distance-calculators.h"
#include "../utils/utils.h" // get_grain

namespace dtwclust {

// =================================================================================================
/* worker to update DTW distance in parallel */
// =================================================================================================

class DtwDistanceUpdater : public RcppParallel::Worker {
public:
    // constructor
    DtwDistanceUpdater(const Rcpp::LogicalVector& id_changed,
                       const Rcpp::IntegerVector& id_nn,
                       Rcpp::NumericMatrix& distmat,
                       const std::shared_ptr<DistanceCalculator>& dist_calculator,
                       int margin)
        : id_changed_(id_changed)
        , id_nn_(id_nn)
        , distmat_(distmat)
        , dist_calculator_(dist_calculator)
        , margin_(margin)
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();
        // update distances
        if (margin_ == 1) {
            for (std::size_t i = begin; i < end; i++) {
                if (id_changed_[i]) {
                    int j = id_nn_[i];
                    distmat_(i,j) = dist_calculator->calculate(i,j);
                }
            }
        }
        else {
            for (std::size_t j = begin; j < end; j++) {
                if (id_changed_[j]) {
                    int i = id_nn_[j];
                    distmat_(i,j) = dist_calculator->calculate(i,j);
                }
            }
        }
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    // input vectors
    const RcppParallel::RVector<int> id_changed_;
    const RcppParallel::RVector<int> id_nn_;
    // output matrix
    RcppParallel::RMatrix<double> distmat_;
    // distance calculator
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    // margin for update
    int margin_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// =================================================================================================
/* find nearest neighbors */
// =================================================================================================

void set_nn(const Rcpp::NumericMatrix& distmat, Rcpp::IntegerVector& nn, const int margin)
{
    if (margin == 1) {
        for (int i = 0; i < distmat.nrow(); i++) {
            double d = distmat(i,0);
            nn[i] = 0;
            for (int j = 1; j < distmat.ncol(); j++) {
                double temp = distmat(i,j);
                if (temp < d) {
                    d = temp;
                    nn[i] = j;
                }
            }
        }
    }
    else {
        for (int j = 0; j < distmat.ncol(); j++) {
            double d = distmat(0,j);
            nn[j] = 0;
            for (int i = 1; i < distmat.nrow(); i++) {
                double temp = distmat(i,j);
                if (temp < d) {
                    d = temp;
                    nn[j] = i;
                }
            }
        }
    }
}

// =================================================================================================
/* check if updates are finished based on indices */
// =================================================================================================

bool check_finished(const Rcpp::IntegerVector& nn,
                    const Rcpp::IntegerVector& nn_prev,
                    Rcpp::LogicalVector& changed)
{
    bool finished = true;
    for (int i = 0; i < nn.length(); i++) {
        if (nn[i] != nn_prev[i]) {
            changed[i] = true;
            finished = false;
        }
        else {
            changed[i] = false;
        }
    }
    return finished;
}

// =================================================================================================
/* main C++ function */
// =================================================================================================

void dtw_lb_cpp(const Rcpp::List& X,
                const Rcpp::List& Y,
                Rcpp::NumericMatrix& distmat,
                const SEXP& DOTS,
                const int margin,
                const int num_threads)
{
    auto dist_calculator = DistanceCalculatorFactory().create("DTW_BASIC", DOTS, X, Y);
    int len = margin == 1 ? distmat.nrow() : distmat.ncol();
    Rcpp::IntegerVector id_nn(len), id_nn_prev(len);
    Rcpp::LogicalVector id_changed(len);
    DtwDistanceUpdater dist_updater(id_changed, id_nn, distmat, dist_calculator, margin);
    set_nn(distmat, id_nn, margin);
    for (int i = 0; i < id_nn.length(); i++) id_nn_prev[i] = id_nn[i] + 1; // initialize different
    int grain = get_grain(len, num_threads);
    while (!check_finished(id_nn, id_nn_prev, id_changed)) {
        Rcpp::checkUserInterrupt();
        // update nn_prev
        for (int i = 0; i < id_nn.length(); i++) id_nn_prev[i] = id_nn[i];
        // calculate dtw distance if necessary
        RcppParallel::parallelFor(0, len, dist_updater, grain);
        // update nearest neighbors
        set_nn(distmat, id_nn, margin);
    }
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

extern "C" SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    Rcpp::NumericMatrix distmat(D);
    dtw_lb_cpp(X, Y, distmat, DOTS, Rcpp::as<int>(MARGIN), Rcpp::as<int>(NUM_THREADS));
    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
