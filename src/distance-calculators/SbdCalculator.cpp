#include "concrete-calculators.h"

#include <RcppArmadillo.h>
#include <RcppParallel.h>

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SbdCalculator::SbdCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
{
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
SbdCalculator* SbdCalculator::clone() const
{
    SbdCalculator* ptr = new SbdCalculator(*this);
    ptr->cc_seq_truncated_ = arma::vec(fftlen_);
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* limits */
// -------------------------------------------------------------------------------------------------

int SbdCalculator::xLimit() const
{
    return x_uv_.length();
}

int SbdCalculator::yLimit() const
{
    return y_uv_.length();
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SbdCalculator::calculate(const int i, const int j)
{
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

} // namespace dtwclust
