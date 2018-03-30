#include "distance-calculators.h"

#include <algorithm> // std::max
#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>

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
DtwBasicCalculator::DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
    , gcm_(nullptr)
{
    Rcpp::List dist_args(DIST_ARGS);
    window_ = Rcpp::as<int>(dist_args["window.size"]);
    norm_ = Rcpp::as<double>(dist_args["norm"]);
    step_ = Rcpp::as<double>(dist_args["step.pattern"]);
    normalize_ = Rcpp::as<bool>(dist_args["normalize"]);
    // set value of max_len_y_
    max_len_y_ = this->maxLength(y_);
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
    return this->calculate(x_[i], y_[j]);
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

double DtwBasicCalculator::calculate(const arma::mat& x, const arma::mat& y) {
    if (!gcm_) return -1;
    return dtw_basic_par(&x[0], &y[0],
                         x.n_rows, y.n_rows, x.n_cols,
                         window_, norm_, step_, normalize_, gcm_);
}

// =================================================================================================
/* Gak */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
GakCalculator::GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
    , logs_(nullptr)
{
    Rcpp::List dist_args(DIST_ARGS);
    sigma_ = Rcpp::as<double>(dist_args["sigma"]);
    window_ = Rcpp::as<int>(dist_args["window.size"]);
    // set values of max_len_*_
    max_len_x_ = this->maxLength(x_);
    max_len_y_ = this->maxLength(y_);
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
    return this->calculate(x_[i], y_[j]);
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

double GakCalculator::calculate(const arma::mat& x, const arma::mat& y) {
    if (!logs_) return -1;
    return logGAK_par(&x[0], &y[0],
                      x.n_rows, y.n_rows, x.n_cols,
                      sigma_, window_, logs_);
}

// =================================================================================================
/* Lbi */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbiCalculator::LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
    , H_(nullptr)
    , L2_(nullptr)
    , U2_(nullptr)
    , LB_(nullptr)
{
    Rcpp::List dist_args(DIST_ARGS);
    p_ = Rcpp::as<int>(dist_args["p"]);
    len_ = Rcpp::as<int>(dist_args["len"]);
    window_ = Rcpp::as<unsigned int>(dist_args["window.size"]);
    Rcpp::List LE((SEXP)dist_args["lower.env"]);
    Rcpp::List UE((SEXP)dist_args["upper.env"]);
    lower_envelopes_ = TSTSList<arma::mat>(LE);
    upper_envelopes_ = TSTSList<arma::mat>(UE);
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
    return this->calculate(x_[i], y_[j], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbiCalculator::calculate(const arma::mat& x, const arma::mat& y,
                                const arma::mat& lower_envelope, const arma::mat& upper_envelope)
{
    return lbi_core(&x[0], &y[0],
                    len_, window_, p_,
                    &lower_envelope[0], &upper_envelope[0],
                    L2_, U2_, H_, LB_);
}

// =================================================================================================
/* Lbk */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbkCalculator::LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , H_(nullptr)
{
    Rcpp::List dist_args(DIST_ARGS);
    p_ = Rcpp::as<int>(dist_args["p"]);
    len_ = Rcpp::as<int>(dist_args["len"]);
    Rcpp::List LE((SEXP)dist_args["lower.env"]);
    Rcpp::List UE((SEXP)dist_args["upper.env"]);
    lower_envelopes_ = TSTSList<arma::mat>(LE);
    upper_envelopes_ = TSTSList<arma::mat>(UE);
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
    return this->calculate(x_[i], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const arma::mat& x,
                                const arma::mat& lower_envelope, const arma::mat& upper_envelope)
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
SbdCalculator::SbdCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
{
    // note cc_seq_truncated_ is not set here, it is allocated for each clone
    Rcpp::List dist_args(DIST_ARGS), x(X), y(Y);
    fftlen_ = Rcpp::as<int>(dist_args["fftlen"]);
    Rcpp::List fftx((SEXP)dist_args["fftx"]);
    Rcpp::List ffty((SEXP)dist_args["ffty"]);
    fftx_ = TSTSList<arma::cx_mat>(fftx);
    ffty_ = TSTSList<arma::cx_mat>(ffty);
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
    return this->calculate(x_[i], y_[j], fftx_[i], ffty_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double SbdCalculator::calculate(const arma::mat& x, const arma::mat& y,
                                const arma::cx_mat& fftx, const arma::cx_mat& ffty)
{
    // already normalizes by length
    arma::mat cc_seq = arma::real(arma::ifft(fftx % ffty));
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
SdtwCalculator::SdtwCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
{
    Rcpp::List dist_args(DIST_ARGS);
    gamma_ = Rcpp::as<double>(dist_args["gamma"]);
    // set values of max_len_*_
    max_len_x_ = this->maxLength(x_);
    max_len_y_ = this->maxLength(y_);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SdtwCalculator::calculate(const int i, const int j) {
    return this->calculate(x_[i], y_[j]);
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

double SdtwCalculator::calculate(const arma::mat& x, const arma::mat& y)
{
    if (!cm_) return -1;
    return sdtw(&x[0], &y[0], x.n_rows, y.n_rows, x.n_cols, gamma_, cm_);
}

} // namespace dtwclust
