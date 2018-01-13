#include "concrete-calculators.h"

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/distances++.h" // lbk_core

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LbkCalculator::LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
    : DistanceCalculator(DIST_ARGS, X, Y)
    , p_(Rcpp::as<int>((SEXP)dist_args_["p"]))
    , len_(Rcpp::as<int>((SEXP)dist_args_["len"]))
{
    x_uv_ = TSTSList<Rcpp::NumericVector>(x_);
    Rcpp::List LE((SEXP)dist_args_["lower.env"]);
    Rcpp::List UE((SEXP)dist_args_["upper.env"]);
    lower_envelopes_ = TSTSList<Rcpp::NumericVector>(LE);
    upper_envelopes_ = TSTSList<Rcpp::NumericVector>(UE);
    H_ = nullptr;
}

// -------------------------------------------------------------------------------------------------
/* destructor */
// -------------------------------------------------------------------------------------------------
LbkCalculator::~LbkCalculator()
{
    if (H_) delete[] H_;
}

// -------------------------------------------------------------------------------------------------
/* clone that sets helper matrix
 *   This is needed because instances of this class are supposed to be called from different
 *   threads, and each one needs its own independent matrix to perform the calculations. Each thread
 *   has to lock a mutex and then call this method before calculating the distance.
 */
// ------------------------------------------------------------------------------------------------
LbkCalculator* LbkCalculator::clone() const
{
    LbkCalculator* ptr = new LbkCalculator(*this);
    ptr->H_ = new double[len_];
    return ptr;
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two lists of series and given indices */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const int i, const int j)
{
    //  is ignored here, only the envelopes matter
    return this->calculate(x_uv_[i], lower_envelopes_[j], upper_envelopes_[j]);
}

// -------------------------------------------------------------------------------------------------
/* compute distance for two series */
// -------------------------------------------------------------------------------------------------
double LbkCalculator::calculate(const RcppParallel::RVector<double>& x,
                                const RcppParallel::RVector<double>& lower_envelope,
                                const RcppParallel::RVector<double>& upper_envelope)
{
    return lbk_core(&x[0], len_, p_, &lower_envelope[0], &upper_envelope[0], H_);
}

} // namespace dtwclust
