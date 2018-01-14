#include "concrete-calculators.h"

#include <algorithm> // std::max

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/distances.h" // logGAK_par

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
GakCalculator::GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
{
    Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
    sigma_ = Rcpp::as<double>(dist_args["sigma"]);
    window_ = Rcpp::as<int>(dist_args["window.size"]);
    is_multivariate_ = Rcpp::as<bool>(dist_args["is.multivariate"]);
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
    // make sure pointer is null
    logs_ = nullptr;
}

// -------------------------------------------------------------------------------------------------
/* destructor */
// -------------------------------------------------------------------------------------------------
GakCalculator::~GakCalculator()
{
    if (logs_) delete[] logs_;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double GakCalculator::calculate(const int i, const int j)
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
GakCalculator* GakCalculator::clone() const
{
    GakCalculator* ptr = new GakCalculator(*this);
    ptr->logs_ = new double[(std::max(max_len_x_, max_len_y_) + 1) * 3];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* limits */
// -------------------------------------------------------------------------------------------------

int GakCalculator::xLimit() const
{
    return is_multivariate_ ? x_mv_.length() : x_uv_.length();
}

int GakCalculator::yLimit() const
{
    return is_multivariate_ ? y_mv_.length() : y_uv_.length();
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------

// univariate
double GakCalculator::calculate(
        const RcppParallel::RVector<double>& x, const RcppParallel::RVector<double>& y)
{
    if (!logs_) return -1;
    int nx = x.length();
    int ny = y.length();
    int num_var = 1;
    return logGAK_par(&x[0], &y[0], nx, ny, num_var, sigma_, window_, logs_);
}

// multivariate
double GakCalculator::calculate(
        const RcppParallel::RMatrix<double>& x, const RcppParallel::RMatrix<double>& y)
{
    if (!logs_) return -1;
    int nx = x.nrow();
    int ny = y.nrow();
    int num_var = x.ncol();
    return logGAK_par(&x[0], &y[0], nx, ny, num_var, sigma_, window_, logs_);
}

} // namespace dtwclust
