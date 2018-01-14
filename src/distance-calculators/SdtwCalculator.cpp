#include "concrete-calculators.h"

#include <utility> // std::move

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/distances++.h" // soft_dtw

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SdtwCalculator::SdtwCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : DistanceCalculator(DIST_ARGS, X, Y)
    , gamma_(Rcpp::as<double>(dist_args_["gamma"]))
    , is_multivariate_(Rcpp::as<bool>(dist_args_["is.multivariate"]))
{
    if (is_multivariate_) {
        x_mv_ = std::move(TSTSList<Rcpp::NumericMatrix>(x_));
        y_mv_ = std::move(TSTSList<Rcpp::NumericMatrix>(y_));
    }
    else {
        x_uv_ = std::move(TSTSList<Rcpp::NumericVector>(x_));
        y_uv_ = std::move(TSTSList<Rcpp::NumericVector>(y_));
    }
    // set values of max_len_*_
    max_len_x_ = this->maxLength(x_, is_multivariate_);
    max_len_y_ = this->maxLength(y_, is_multivariate_);
    // make sure pointer is null
    cm_ = nullptr;
}

// -------------------------------------------------------------------------------------------------
/* destructor */
// -------------------------------------------------------------------------------------------------
SdtwCalculator::~SdtwCalculator()
{
    if (cm_) delete[] cm_;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SdtwCalculator::calculate(const int i, const int j)
{
    if (is_multivariate_) {
        return this->calculate(x_mv_[i], y_mv_[j]);
    }
    else {
        return this->calculate(x_uv_[i], y_uv_[j]);
    }
}

// -------------------------------------------------------------------------------------------------
/* clone that sets helper matrix
 *   This is needed because instances of this class are supposed to be called from different
 *   threads, and each one needs its own independent matrix to perform the calculations. Each thread
 *   has to lock a mutex and then call this method before calculating the distance.
 */
// ------------------------------------------------------------------------------------------------
SdtwCalculator* SdtwCalculator::clone() const
{
    SdtwCalculator* ptr = new SdtwCalculator(*this);
    ptr->cm_ = new double[(max_len_x_ + 1) * (max_len_y_ + 1)];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------

// univariate
double SdtwCalculator::calculate(
        const RcppParallel::RVector<double>& x, const RcppParallel::RVector<double>& y)
{
    if (!cm_) return -1;
    int nx = x.length();
    int ny = y.length();
    int num_var = 1;
    return sdtw(&x[0], &y[0], nx, ny, num_var, gamma_, cm_);
}

// multivariate
double SdtwCalculator::calculate(
        const RcppParallel::RMatrix<double>& x, const RcppParallel::RMatrix<double>& y)
{
    if (!cm_) return -1;
    int nx = x.nrow();
    int ny = y.nrow();
    int num_var = x.ncol();
    return sdtw(&x[0], &y[0], nx, ny, num_var, gamma_, cm_);
}

} // namespace dtwclust
