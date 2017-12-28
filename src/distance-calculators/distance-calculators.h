#ifndef _DTWCLUST_DISTANCE_CALCULATORS_HPP_
#define _DTWCLUST_DISTANCE_CALCULATORS_HPP_

#include <RcppArmadillo.h> // arma:: referenced here
#include <memory> // *_ptr
#include <string>

namespace dtwclust {

// =================================================================================================
/* DistanceCalculator (base + factory + concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distance calculator */
// -------------------------------------------------------------------------------------------------
class DistanceCalculator
{
public:
    virtual ~DistanceCalculator() {};
    virtual double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) = 0;
protected:
    DistanceCalculator(const SEXP& DIST_ARGS) : dist_args_(DIST_ARGS) {};
    Rcpp::List dist_args_;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistanceCalculatorFactory
{
public:
    std::shared_ptr<DistanceCalculator> create(const SEXP& DIST, const SEXP& DIST_ARGS);
    std::shared_ptr<DistanceCalculator> create(const std::string& dist, const SEXP& DIST_ARGS);
};

// -------------------------------------------------------------------------------------------------
/* dtw_basic calculator */
// -------------------------------------------------------------------------------------------------
class DtwBasicDistanceCalculator : public DistanceCalculator
{
public:
    DtwBasicDistanceCalculator(const SEXP& DIST_ARGS);
    double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                     const int i, const int j) override;
private:
    double calculate(const SEXP& X, const SEXP& Y);
    SEXP window_, norm_, step_, backtrack_, gcm_;
    bool is_multivariate_, normalize_;
};

// -------------------------------------------------------------------------------------------------
/* lb_keogh calculator */
// -------------------------------------------------------------------------------------------------
class LbkDistanceCalculator : public DistanceCalculator
{
public:
    LbkDistanceCalculator(const SEXP& DIST_ARGS);
    double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                     const int i, const int j) override;
private:
    double calculate(const Rcpp::NumericVector& x,
                     const Rcpp::NumericVector& lower_envelope,
                     const Rcpp::NumericVector& upper_envelope);
    Rcpp::List lower_envelopes_, upper_envelopes_;
    Rcpp::NumericVector H_;
    int p_;
};

// -------------------------------------------------------------------------------------------------
/* lb_improved calculator */
// -------------------------------------------------------------------------------------------------
class LbiDistanceCalculator : public DistanceCalculator
{
public:
    LbiDistanceCalculator(const SEXP& DIST_ARGS);
    double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                     const int i, const int j) override;
private:
    double calculate(const Rcpp::NumericVector& x,
                     const Rcpp::NumericVector& y,
                     const Rcpp::NumericVector& lower_envelope,
                     const Rcpp::NumericVector& upper_envelope);
    Rcpp::List lower_envelopes_, upper_envelopes_;
    Rcpp::NumericVector H_, L2_, U2_, LB_;
    unsigned int window_size_;
    int p_;
};

// -------------------------------------------------------------------------------------------------
/* soft-DTW calculator */
// -------------------------------------------------------------------------------------------------
class SdtwDistanceCalculator : public DistanceCalculator
{
public:
    SdtwDistanceCalculator(const SEXP& DIST_ARGS);
    double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                     const int i, const int j) override;
private:
    double calculate(const SEXP& X, const SEXP& Y);
    SEXP gamma_, costmat_, mv_;
};

// -------------------------------------------------------------------------------------------------
/* gak calculator */
// -------------------------------------------------------------------------------------------------
class GakDistanceCalculator : public DistanceCalculator
{
public:
    GakDistanceCalculator(const SEXP& DIST_ARGS);
    double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                     const int i, const int j) override;
private:
    double calculate(const SEXP& X, const SEXP& Y);
    SEXP sigma_, window_, logs_;
    bool is_multivariate_;
};

// -------------------------------------------------------------------------------------------------
/* sbd calculator */
// -------------------------------------------------------------------------------------------------
class SbdDistanceCalculator : public DistanceCalculator
{
public:
    SbdDistanceCalculator(const SEXP& DIST_ARGS);
    double calculate(const Rcpp::List& X, const Rcpp::List& Y,
                     const int i, const int j) override;
private:
    double calculate(const arma::vec& x, const arma::vec& y,
                     const arma::cx_vec& fftx, const arma::cx_vec& ffty);
    arma::vec cc_seq_truncated_;
    Rcpp::List fftx_, ffty_;
    int fftlen_;
};

} // namespace dtwclust

#endif // _DTWCLUST_DISTANCE_CALCULATORS_HPP_
