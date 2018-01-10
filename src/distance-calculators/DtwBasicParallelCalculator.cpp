#include "distance-calculators.h"

#include <vector>

#include <RcppParallel.h>

#include "../distances/distances.h" // dtw_basic_par

namespace dtwclust {

// =================================================================================================
/* dtw_basic distance parallel calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
DtwBasicParallelCalculator::DtwBasicParallelCalculator(
    const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : DistanceCalculator(DIST_ARGS, X, Y)
    , window_(Rcpp::as<int>(dist_args_["window.size"]))
    , norm_(Rcpp::as<double>(dist_args_["norm"]))
    , step_(Rcpp::as<double>(dist_args_["step.pattern"]))
    , is_multivariate_(Rcpp::as<bool>(dist_args_["is.multivariate"]))
{
    if (is_multivariate_) { // nocov start
        x_mv_ = TSTSList<Rcpp::NumericMatrix>(x_);
        y_mv_ = TSTSList<Rcpp::NumericMatrix>(y_);
    } // nocov end
    else {
        x_uv_ = TSTSList<Rcpp::NumericVector>(x_);
        y_uv_ = TSTSList<Rcpp::NumericVector>(y_);
    }
}

// -------------------------------------------------------------------------------------------------
/* set pointer to helper matrix */
// -------------------------------------------------------------------------------------------------
void DtwBasicParallelCalculator::setGcm(double * const gcm)
{ gcm_ = gcm; }

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------

// univariate
double DtwBasicParallelCalculator::calculate(
        const RcppParallel::RVector<double>& x, const RcppParallel::RVector<double>& y)
{
    if (gcm_ == nullptr) return -1;
    int nx = x.length();
    int ny = y.length();
    int num_var = 1;
    return dtw_basic_par(&x[0], &y[0], nx, ny, num_var, window_, norm_, step_, gcm_);
}

// multivariate
double DtwBasicParallelCalculator::calculate( // nocov start
        const RcppParallel::RMatrix<double>& x, const RcppParallel::RMatrix<double>& y)
{
    if (gcm_ == nullptr) return -1;
    int nx = x.nrow();
    int ny = y.nrow();
    int num_var = x.ncol();
    return dtw_basic_par(&x[0], &y[0], nx, ny, num_var, window_, norm_, step_, gcm_);
} // nocov end

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double DtwBasicParallelCalculator::calculate(const int i, const int j)
{
    if (is_multivariate_) {
        return this->calculate(x_mv_[i], y_mv_[j]); // nocov
    }
    else {
        return this->calculate(x_uv_[i], y_uv_[j]);
    }
}

} // namespace dtwclust
