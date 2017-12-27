#include "dtwclust++.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* sbd distance calculator */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SbdDistanceCalculator::SbdDistanceCalculator(const SEXP& DIST_ARGS)
    : DistanceCalculator(DIST_ARGS)
{
    fftlen_ = Rcpp::as<int>((SEXP)dist_args_["fftlen"]);
    fftx_ = dist_args_["fftx"];
    ffty_ = dist_args_["ffty"];
    cc_seq_truncated_ = arma::vec(fftlen_);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double SbdDistanceCalculator::calculateDistance(const arma::vec& x, const arma::vec& y,
                                                const arma::cx_vec& fftx, const arma::cx_vec& ffty)
{
    // already normalizes by length
    arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));
    double x_norm = arma::norm(x);
    double y_norm = arma::norm(y);

    // reorder truncated sequence
    int id = 0;
    for (unsigned int i = fftlen_ - y.size() + 1; i < cc_seq.size(); i++) {
        cc_seq_truncated_[id] = cc_seq[i];
        id++;
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        cc_seq_truncated_[id] = cc_seq[i];
        id++;
    }

    // get max
    double cc_max = R_NegInf;
    double den = x_norm * y_norm;
    for (int i = 0; i < id; i++) {
        double this_cc = cc_seq_truncated_[i] / den;
        if (this_cc > cc_max) cc_max = this_cc;
    }

    return 1 - cc_max;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double SbdDistanceCalculator::calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                                                const int i, const int j)
{
    // in two steps to avoid disambiguation
    Rcpp::NumericVector x_rcpp(X[i]);
    Rcpp::NumericVector y_rcpp(Y[j]);
    Rcpp::ComplexVector fftx_rcpp(fftx_[i]);
    Rcpp::ComplexVector ffty_rcpp(ffty_[j]);
    arma::vec x(x_rcpp), y(y_rcpp);
    arma::cx_vec fftx(fftx_rcpp), ffty(ffty_rcpp);
    return this->calculateDistance(x, y, fftx, ffty);
}

} // namespace dtwclust
