#include "concrete-calculators.h"

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/distances.h" // dtw_basic_par

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
DtwBasicCalculator::DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
{
    Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
    window_ = Rcpp::as<int>(dist_args["window.size"]);
    norm_ = Rcpp::as<double>(dist_args["norm"]);
    step_ = Rcpp::as<double>(dist_args["step.pattern"]);
    normalize_ = Rcpp::as<bool>(dist_args["normalize"]);
    is_multivariate_ = Rcpp::as<bool>(dist_args["is.multivariate"]);
    if (is_multivariate_) {
        x_mv_ = TSTSList<Rcpp::NumericMatrix>(x);
        y_mv_ = TSTSList<Rcpp::NumericMatrix>(y);
    }
    else {
        x_uv_ = TSTSList<Rcpp::NumericVector>(x);
        y_uv_ = TSTSList<Rcpp::NumericVector>(y);
    }
    // set value of max_len_y_
    max_len_y_ = this->maxLength(y, is_multivariate_);
    // make sure pointer is null
    gcm_ = nullptr;
}

// -------------------------------------------------------------------------------------------------
/* destructor */
// -------------------------------------------------------------------------------------------------
DtwBasicCalculator::~DtwBasicCalculator()
{
    if (gcm_) delete[] gcm_;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double DtwBasicCalculator::calculate(const int i, const int j)
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
DtwBasicCalculator* DtwBasicCalculator::clone() const
{
    DtwBasicCalculator* ptr = new DtwBasicCalculator(*this);
    ptr->gcm_ = new double[2 * (max_len_y_ + 1)];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* limits */
// -------------------------------------------------------------------------------------------------

int DtwBasicCalculator::xLimit() const
{
    return is_multivariate_ ? x_mv_.length() : x_uv_.length();
}

int DtwBasicCalculator::yLimit() const
{
    return is_multivariate_ ? y_mv_.length() : y_uv_.length();
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------

// univariate
double DtwBasicCalculator::calculate(
        const RcppParallel::RVector<double>& x, const RcppParallel::RVector<double>& y)
{
    if (!gcm_) return -1;
    int nx = x.length();
    int ny = y.length();
    int num_var = 1;
    return dtw_basic_par(&x[0], &y[0],
                         nx, ny, num_var,
                         window_, norm_, step_, normalize_,
                         gcm_, false, nullptr, nullptr, nullptr);
}

// multivariate
double DtwBasicCalculator::calculate(
        const RcppParallel::RMatrix<double>& x, const RcppParallel::RMatrix<double>& y)
{
    if (!gcm_) return -1;
    int nx = x.nrow();
    int ny = y.nrow();
    int num_var = x.ncol();
    return dtw_basic_par(&x[0], &y[0],
                         nx, ny, num_var,
                         window_, norm_, step_, normalize_,
                         gcm_, false, nullptr, nullptr, nullptr);
}

} // namespace dtwclust
