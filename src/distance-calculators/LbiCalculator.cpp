#include "concrete-calculators.h"

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/distances.h" // lbi_core

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbiCalculator::LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
{
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
LbiCalculator::~LbiCalculator()
{
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
LbiCalculator* LbiCalculator::clone() const
{
    LbiCalculator* ptr = new LbiCalculator(*this);
    ptr->H_ = new double[len_];
    ptr->L2_ = new double[len_];
    ptr->U2_ = new double[len_];
    ptr->LB_ = new double[len_];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* limits */
// -------------------------------------------------------------------------------------------------

int LbiCalculator::xLimit() const
{
    return x_uv_.length();
}

int LbiCalculator::yLimit() const
{
    return y_uv_.length();
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double LbiCalculator::calculate(const int i, const int j)
{
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

} // namespace dtwclust
