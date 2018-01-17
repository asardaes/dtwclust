#ifndef DTWCLUST_CONCRETE_CALCULATORS_HPP_
#define DTWCLUST_CONCRETE_CALCULATORS_HPP_

#include <RcppArmadillo.h> // arma:: referenced here
#include <RcppParallel.h>

#include "../utils/TSTSList.h"
#include "distance-calculators.h"

namespace dtwclust {

// =================================================================================================
/* DistanceCalculator (concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* dtw_basic calculator */
// -------------------------------------------------------------------------------------------------
class DtwBasicCalculator : public DistanceCalculator
{
public:
    DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    ~DtwBasicCalculator();
    double calculate(const int i, const int j) override;
    DtwBasicCalculator* clone() const override;
    int xLimit() const override;
    int yLimit() const override;

private:
    // method calculate
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& y);
    double calculate(const RcppParallel::RMatrix<double>& x,
                     const RcppParallel::RMatrix<double>& y);
    // input parameters
    int window_;
    double norm_, step_;
    bool normalize_, is_multivariate_;
    // input series (univariate)
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    // input series (multivariate)
    TSTSList<Rcpp::NumericMatrix> x_mv_, y_mv_;
    // helper "matrix"
    double* gcm_;
    // to dimension gcm_
    int max_len_y_;
};

// -------------------------------------------------------------------------------------------------
/* gak calculator */
// -------------------------------------------------------------------------------------------------
class GakCalculator : public DistanceCalculator
{
public:
    GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    ~GakCalculator();
    double calculate(const int i, const int j) override;
    GakCalculator* clone() const override;
    int xLimit() const override;
    int yLimit() const override;

private:
    // method calculate
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& y);
    double calculate(const RcppParallel::RMatrix<double>& x,
                     const RcppParallel::RMatrix<double>& y);
    // input parameters
    double sigma_;
    int window_;
    bool is_multivariate_;
    // input series (univariate)
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    // input series (multivariate)
    TSTSList<Rcpp::NumericMatrix> x_mv_, y_mv_;
    // helper "matrix"
    double* logs_;
    // to dimension logs_
    int max_len_x_, max_len_y_;
};

// -------------------------------------------------------------------------------------------------
/* lb_improved calculator */
// -------------------------------------------------------------------------------------------------
class LbiCalculator : public DistanceCalculator
{
public:
    LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    ~LbiCalculator();
    double calculate(const int i, const int j) override;
    LbiCalculator* clone() const override;
    int xLimit() const override;
    int yLimit() const override;

private:
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& y,
                     const RcppParallel::RVector<double>& lower_envelope,
                     const RcppParallel::RVector<double>& upper_envelope);
    int p_, len_;
    unsigned int window_;
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_, lower_envelopes_, upper_envelopes_;
    double *H_, *L2_, *U2_, *LB_;
};

// -------------------------------------------------------------------------------------------------
/* lb_keogh calculator */
// -------------------------------------------------------------------------------------------------
class LbkCalculator : public DistanceCalculator
{
public:
    LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    ~LbkCalculator();
    double calculate(const int i, const int j) override;
    LbkCalculator* clone() const override;
    int xLimit() const override;
    int yLimit() const override;

private:
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& lower_envelope,
                     const RcppParallel::RVector<double>& upper_envelope);
    int p_, len_;
    TSTSList<Rcpp::NumericVector> x_uv_, lower_envelopes_, upper_envelopes_;
    double* H_;
};

// -------------------------------------------------------------------------------------------------
/* sbd calculator */
// -------------------------------------------------------------------------------------------------
class SbdCalculator : public DistanceCalculator
{
public:
    SbdCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    SbdCalculator* clone() const override;
    int xLimit() const override;
    int yLimit() const override;

private:
    double calculate(const arma::vec& x, const arma::vec& y,
                     const arma::cx_vec& fftx, const arma::cx_vec& ffty);
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    TSTSList<Rcpp::ComplexVector> fftx_, ffty_;
    arma::vec cc_seq_truncated_;
    int fftlen_;
};

// -------------------------------------------------------------------------------------------------
/* soft-DTW calculator */
// -------------------------------------------------------------------------------------------------
class SdtwCalculator : public DistanceCalculator
{
public:
    SdtwCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    ~SdtwCalculator();
    double calculate(const int i, const int j) override;
    SdtwCalculator* clone() const override;
    int xLimit() const override;
    int yLimit() const override;

private:
    // method calculate
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& y);
    double calculate(const RcppParallel::RMatrix<double>& x,
                     const RcppParallel::RMatrix<double>& y);
    // input parameters
    double gamma_;
    bool is_multivariate_;
    // input series (univariate)
    TSTSList<Rcpp::NumericVector> x_uv_, y_uv_;
    // input series (multivariate)
    TSTSList<Rcpp::NumericMatrix> x_mv_, y_mv_;
    // helper "matrix"
    double* cm_;
    // to dimension cm_
    int max_len_x_, max_len_y_;
};

} // namespace dtwclust

#endif // DTWCLUST_CONCRETE_CALCULATORS_HPP_
