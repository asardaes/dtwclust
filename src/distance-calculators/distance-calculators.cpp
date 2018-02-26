#include "distance-calculators.h"

#include <algorithm> // std::max
#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/distances-details.h"  // dtw_basic_par, logGAK_par, lbi_core, lbk_core, soft_dtw
#include "../utils/SurrogateMatrix.h"

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::create(const SEXP& DIST, const SEXP& DIST_ARGS,
                                  const SEXP& X, const SEXP& Y)
{
    std::string dist = Rcpp::as<std::string>(DIST);
    return this->create(dist, DIST_ARGS, X, Y);
}

std::shared_ptr<DistanceCalculator>
DistanceCalculatorFactory::create(const std::string& dist, const SEXP& DIST_ARGS,
                                  const SEXP& X, const SEXP& Y)
{
    if (dist == "DTW_BASIC")
        return std::make_shared<DtwBasicCalculator>(DIST_ARGS, X, Y);
    else if (dist == "LBK")
        return std::make_shared<LbkCalculator>(DIST_ARGS, X, Y);
    else if (dist == "LBI")
        return std::make_shared<LbiCalculator>(DIST_ARGS, X, Y);
    else if (dist == "SDTW")
        return std::make_shared<SdtwCalculator>(DIST_ARGS, X, Y);
    else if (dist == "GAK")
        return std::make_shared<GakCalculator>(DIST_ARGS, X, Y);
    else if (dist == "SBD")
        return std::make_shared<SbdCalculator>(DIST_ARGS, X, Y);
    else
        Rcpp::stop("Unknown distance measure"); // nocov
}

// =================================================================================================
/* DtwBasic */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
DtwBasicCalculator::DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y) {
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
DtwBasicCalculator::~DtwBasicCalculator() {
    if (gcm_) delete[] gcm_;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double DtwBasicCalculator::calculate(const int i, const int j) {
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
DtwBasicCalculator* DtwBasicCalculator::clone() const {
    DtwBasicCalculator* ptr = new DtwBasicCalculator(*this);
    ptr->gcm_ = new double[2 * (max_len_y_ + 1)];
    return ptr;
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

// =================================================================================================
/* Gak */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
GakCalculator::GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y) {
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
GakCalculator::~GakCalculator() {
    if (logs_) delete[] logs_;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double GakCalculator::calculate(const int i, const int j) {
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
GakCalculator* GakCalculator::clone() const {
    GakCalculator* ptr = new GakCalculator(*this);
    ptr->logs_ = new double[(std::max(max_len_x_, max_len_y_) + 1) * 3];
    return ptr;
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

// =================================================================================================
/* Lbi */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbiCalculator::LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y) {
    Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
    p_ = Rcpp::as<int>(dist_args["p"]);
    len_ = Rcpp::as<int>(dist_args["len"]);
    window_ = Rcpp::as<unsigned int>(dist_args["window.size"]);
    x_uv_ = TSTSList<Rcpp::NumericVector>(x);
    y_uv_ = TSTSList<Rcpp::NumericVector>(y);
    Rcpp::List LE((SEXP)dist_args["lower.env"]);
    Rcpp::List UE((SEXP)dist_args["upper.env"]);
    lower_envelopes_ = TSTSList<Rcpp::NumericVector>(LE);
    upper_envelopes_ = TSTSList<Rcpp::NumericVector>(UE);
    H_ = nullptr;
    L2_ = nullptr;
    U2_ = nullptr;
    LB_ = nullptr;
}

// -------------------------------------------------------------------------------------------------
/* destructor */
// -------------------------------------------------------------------------------------------------
LbiCalculator::~LbiCalculator() {
    if (H_) delete[] H_;
    if (L2_) delete[] L2_;
    if (U2_) delete[] U2_;
    if (LB_) delete[] LB_;
}

// -------------------------------------------------------------------------------------------------
/* clone that sets helper matrix
*   This is needed because instances of this class are supposed to be called from different
*   threads, and each one needs its own independent matrix to perform the calculations. Each thread
*   has to lock a mutex and then call this method before calculating the distance.
*/
// ------------------------------------------------------------------------------------------------
LbiCalculator* LbiCalculator::clone() const {
    LbiCalculator* ptr = new LbiCalculator(*this);
    ptr->H_ = new double[len_];
    ptr->L2_ = new double[len_];
    ptr->U2_ = new double[len_];
    ptr->LB_ = new double[len_];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double LbiCalculator::calculate(const int i, const int j) {
    return this->calculate(x_uv_[i], y_uv_[j], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbiCalculator::calculate(const RcppParallel::RVector<double>& x,
                                const RcppParallel::RVector<double>& y,
                                const RcppParallel::RVector<double>& lower_envelope,
                                const RcppParallel::RVector<double>& upper_envelope)
{
    return lbi_core(&x[0], &y[0], len_, window_, p_, &lower_envelope[0], &upper_envelope[0],
                    L2_, U2_, H_, LB_);
}

// =================================================================================================
/* Lbk */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbkCalculator::LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y) {
    Rcpp::List dist_args(DIST_ARGS), x(X);
    p_ = Rcpp::as<int>(dist_args["p"]);
    len_ = Rcpp::as<int>(dist_args["len"]);
    x_uv_ = TSTSList<Rcpp::NumericVector>(x);
    Rcpp::List LE((SEXP)dist_args["lower.env"]);
    Rcpp::List UE((SEXP)dist_args["upper.env"]);
    lower_envelopes_ = TSTSList<Rcpp::NumericVector>(LE);
    upper_envelopes_ = TSTSList<Rcpp::NumericVector>(UE);
    H_ = nullptr;
}

// -------------------------------------------------------------------------------------------------
/* destructor */
// -------------------------------------------------------------------------------------------------
LbkCalculator::~LbkCalculator() {
    if (H_) delete[] H_;
}

// -------------------------------------------------------------------------------------------------
/* clone that sets helper matrix
*   This is needed because instances of this class are supposed to be called from different
*   threads, and each one needs its own independent matrix to perform the calculations. Each thread
*   has to lock a mutex and then call this method before calculating the distance.
*/
// ------------------------------------------------------------------------------------------------
LbkCalculator* LbkCalculator::clone() const {
    LbkCalculator* ptr = new LbkCalculator(*this);
    ptr->H_ = new double[len_];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const int i, const int j) {
    // y is ignored here, only the envelopes matter
    return this->calculate(x_uv_[i], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const RcppParallel::RVector<double>& x,
                                const RcppParallel::RVector<double>& lower_envelope,
                                const RcppParallel::RVector<double>& upper_envelope)
{
    if (!H_) return -1;
    return lbk_core(&x[0], len_, p_, &lower_envelope[0], &upper_envelope[0], H_);
}

// =================================================================================================
/* Sbd */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SbdCalculator::SbdCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y) {
    // note cc_seq_truncated_ is not set here, it is allocated for each clone
    Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
    fftlen_ = Rcpp::as<int>(dist_args["fftlen"]);
    x_uv_ = TSTSList<Rcpp::NumericVector>(x);
    y_uv_ = TSTSList<Rcpp::NumericVector>(y);
    Rcpp::List fftx((SEXP)dist_args["fftx"]);
    Rcpp::List ffty((SEXP)dist_args["ffty"]);
    fftx_ = TSTSList<Rcpp::ComplexVector>(fftx);
    ffty_ = TSTSList<Rcpp::ComplexVector>(ffty);
}

// -------------------------------------------------------------------------------------------------
/* clone */
// ------------------------------------------------------------------------------------------------
SbdCalculator* SbdCalculator::clone() const {
    SbdCalculator* ptr = new SbdCalculator(*this);
    ptr->cc_seq_truncated_ = arma::vec(fftlen_);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SbdCalculator::calculate(const int i, const int j) {
    RcppParallel::RVector<double>& x_rcpp = x_uv_[i];
    RcppParallel::RVector<double>& y_rcpp = y_uv_[j];
    // avoid copying data (ptr, length, copy, strict)
    const arma::vec x(x_rcpp.begin(), x_rcpp.length(), false, true);
    const arma::vec y(y_rcpp.begin(), y_rcpp.length(), false, true);
    // calculate distance
    return this->calculate(x, y, fftx_[i], ffty_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double SbdCalculator::calculate(const arma::vec& x, const arma::vec& y,
                                const arma::cx_vec& fftx, const arma::cx_vec& ffty)
{
    // already normalizes by length
    arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));
    double x_norm = arma::norm(x);
    double y_norm = arma::norm(y);
    // reorder truncated sequence
    int id = 0;
    for (unsigned int i = fftlen_ - y.size() + 1; i < cc_seq.size(); i++) {
        cc_seq_truncated_[id++] = cc_seq[i];
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        cc_seq_truncated_[id++] = cc_seq[i];
    }
    // get max
    double cc_max = R_NegInf;
    double den = x_norm * y_norm;
    for (int i = 0; i < id; i++) {
        double this_cc = cc_seq_truncated_[i] / den;
        if (this_cc > cc_max) cc_max = this_cc;
    }
    // return distance value
    return 1 - cc_max;
}

// =================================================================================================
/* Sdtw */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SdtwCalculator::SdtwCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y) {
    Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
    gamma_ = Rcpp::as<double>(dist_args["gamma"]);
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
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SdtwCalculator::calculate(const int i, const int j) {
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
SdtwCalculator* SdtwCalculator::clone() const {
    SdtwCalculator* ptr = new SdtwCalculator(*this);
    ptr->cm_ = SurrogateMatrix<double>(max_len_x_ + 2, max_len_y_ + 2);
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
    return sdtw(&x[0], &y[0], nx, ny, 1, gamma_, cm_);
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
