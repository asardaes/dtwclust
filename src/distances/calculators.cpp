#include "calculators.h"

#include <algorithm> // std::max
#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>

#include "details.h"  // dtw_basic, logGAK, lbi_core, lbk_core, soft_dtw
#include "../utils/SurrogateMatrix.h"
#include "../utils/utils.h" // id_t

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
DtwBasicCalculator::DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
{
    Rcpp::List dist_args(DIST_ARGS);
    window_ = Rcpp::as<int>(dist_args["window.size"]);
    norm_ = Rcpp::as<double>(dist_args["norm"]);
    step_ = Rcpp::as<double>(dist_args["step.pattern"]);
    normalize_ = Rcpp::as<bool>(dist_args["normalize"]);
    sqrt_dist_ = Rcpp::as<bool>(dist_args["sqrt.dist"]);
    // set value of max_len_y_
    max_len_y_ = this->maxLength(y_);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
double DtwBasicCalculator::calculate(const id_t i, const id_t j) {
    return this->calculate(x_[i], y_[j]);
}

// -------------------------------------------------------------------------------------------------
/* clone */
DtwBasicCalculator* DtwBasicCalculator::clone() const {
    DtwBasicCalculator* ptr = new DtwBasicCalculator(*this);
    ptr->lcm_ = SurrogateMatrix<double>(2, max_len_y_ + 1);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
double DtwBasicCalculator::calculate(const arma::mat& x, const arma::mat& y) {
    if (!lcm_) return -1;

    SurrogateMatrix<const double> temp_x(x.n_rows, x.n_cols, &x[0]);
    SurrogateMatrix<const double> temp_y(y.n_rows, y.n_cols, &y[0]);
    return dtw_basic(lcm_, temp_x, temp_y, window_, norm_, step_, normalize_, sqrt_dist_);
}

// =================================================================================================
/* Gak */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
GakCalculator::GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
{
    Rcpp::List dist_args(DIST_ARGS);
    sigma_ = Rcpp::as<double>(dist_args["sigma"]);
    window_ = Rcpp::as<int>(dist_args["window.size"]);
    // set values of max_len_*_
    max_len_x_ = this->maxLength(x_);
    max_len_y_ = this->maxLength(y_);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
double GakCalculator::calculate(const id_t i, const id_t j) {
    return this->calculate(x_[i], y_[j]);
}

// -------------------------------------------------------------------------------------------------
/* clone */
GakCalculator* GakCalculator::clone() const {
    GakCalculator* ptr = new GakCalculator(*this);
    ptr->logs_ = SurrogateMatrix<double>(std::max(max_len_x_, max_len_y_) + 1, 3);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
double GakCalculator::calculate(const arma::mat& x, const arma::mat& y) {
    if (!logs_) return -1;

    SurrogateMatrix<const double> temp_x(x.n_rows, x.n_cols, &x[0]);
    SurrogateMatrix<const double> temp_y(y.n_rows, y.n_cols, &y[0]);
    return logGAK_c(temp_x, temp_y, sigma_, static_cast<id_t>(window_), logs_);
}

// =================================================================================================
/* Lbi */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
LbiCalculator::LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
    , y_(Y)
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
/* clone */
LbiCalculator* LbiCalculator::clone() const {
    LbiCalculator* ptr = new LbiCalculator(*this);
    ptr->H_ = SurrogateMatrix<double>(len_, 1);
    ptr->L2_ = SurrogateMatrix<double>(len_, 1);
    ptr->U2_ = SurrogateMatrix<double>(len_, 1);
    ptr->LB_ = SurrogateMatrix<double>(len_, 1);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
double LbiCalculator::calculate(const id_t i, const id_t j) {
    return this->calculate(x_[i], y_[j], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
double LbiCalculator::calculate(const arma::mat& x, const arma::mat& y,
                                const arma::mat& lower_envelope, const arma::mat& upper_envelope)
{
    SurrogateMatrix<const double> temp_x(len_, 1, &x[0]);
    SurrogateMatrix<const double> temp_y(len_, 1, &y[0]);
    SurrogateMatrix<const double> temp_l(len_, 1, &lower_envelope[0]);
    SurrogateMatrix<const double> temp_u(len_, 1, &upper_envelope[0]);
    return lbi_core(temp_x, temp_y,
                    window_, p_,
                    temp_l, temp_u,
                    L2_, U2_, H_, LB_);
}

// =================================================================================================
/* Lbk */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
LbkCalculator::LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : x_(X)
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
/* clone */
LbkCalculator* LbkCalculator::clone() const {
    LbkCalculator* ptr = new LbkCalculator(*this);
    ptr->H_ = SurrogateMatrix<double>(len_, 1);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
double LbkCalculator::calculate(const id_t i, const id_t j) {
    // y is ignored here, only the envelopes matter
    return this->calculate(x_[i], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
double LbkCalculator::calculate(const arma::mat& x,
                                const arma::mat& lower_envelope, const arma::mat& upper_envelope)
{
    if (!H_) return -1;

    SurrogateMatrix<const double> temp_x(len_, 1, &x[0]);
    SurrogateMatrix<const double> temp_l(len_, 1, &lower_envelope[0]);
    SurrogateMatrix<const double> temp_u(len_, 1, &upper_envelope[0]);
    return lbk_core(temp_x, p_, temp_l, temp_u, H_);
}

// =================================================================================================
/* Sbd */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
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
SbdCalculator* SbdCalculator::clone() const {
    SbdCalculator* ptr = new SbdCalculator(*this);
    ptr->cc_seq_truncated_ = arma::vec(fftlen_);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
double SbdCalculator::calculate(const id_t i, const id_t j) {
    return this->calculate(x_[i], y_[j], fftx_[i], ffty_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
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
double SdtwCalculator::calculate(const id_t i, const id_t j) {
    return this->calculate(x_[i], y_[j]);
}

// -------------------------------------------------------------------------------------------------
/* clone */
SdtwCalculator* SdtwCalculator::clone() const {
    SdtwCalculator* ptr = new SdtwCalculator(*this);
    ptr->cm_ = SurrogateMatrix<double>(max_len_x_ + 2, max_len_y_ + 2);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
double SdtwCalculator::calculate(const arma::mat& x, const arma::mat& y)
{
    if (!cm_) return -1;

    SurrogateMatrix<const double> temp_x(x.n_rows, x.n_cols, &x[0]);
    SurrogateMatrix<const double> temp_y(y.n_rows, y.n_cols, &y[0]);
    return sdtw(temp_x, temp_y, gamma_, cm_);
}

} // namespace dtwclust
