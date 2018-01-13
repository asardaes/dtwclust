#ifndef DTWCLUST_DISTANCE_CALCULATORS_HPP_
#define DTWCLUST_DISTANCE_CALCULATORS_HPP_

#include <complex>
#include <memory> // *_ptr
#include <string>
#include <type_traits> // conditional, is_same
#include <vector>

#include <RcppArmadillo.h> // arma:: referenced here
#include <RcppParallel.h>

namespace dtwclust {

// =================================================================================================
/* Thread-Safe Time-Series List
 *   It should be constructed outside the thread though.
 */
// =================================================================================================

// primary one
template<typename SeriesType>
class TSTSList
{
public:
    typedef typename std::conditional<
        std::is_same<SeriesType, Rcpp::NumericVector>::value,
        RcppParallel::RVector<double>,
        RcppParallel::RMatrix<double>
    >::type TSTSType;
    // constructors
    TSTSList() {} // dummy
    TSTSList(const Rcpp::List& series)
    {
        for (const SEXP& x : series) {
            SeriesType x_rcpp(x);
            series_.push_back(TSTSType(x_rcpp));
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return series_[i]; }
    const TSTSType& operator[](const int i) const { return series_[i]; }

private:
    std::vector<TSTSType> series_;
};

// specialization for complex vectors
template<>
class TSTSList<Rcpp::ComplexVector>
{
public:
    typedef arma::cx_vec TSTSType;
    // constructors
    TSTSList() {} // dummy
    TSTSList(const Rcpp::List& series)
    {
        for (const SEXP& x : series) {
            Rcpp::ComplexVector x_rcpp(x);
            // see http://rcpp-devel.r-forge.r-project.narkive.com/o5ubHVos/multiplication-of-complexvector
            arma::cx_vec x_arma(reinterpret_cast<std::complex<double>*>(x_rcpp.begin()),
                                x_rcpp.length(), false, true);
            series_.push_back(x_arma);
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return series_[i]; }
    const TSTSType& operator[](const int i) const { return series_[i]; }

private:
    std::vector<TSTSType> series_;
};

// =================================================================================================
/* DistanceCalculator (base + factory + concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distance calculator */
// -------------------------------------------------------------------------------------------------
class DistanceCalculator
{
public:
    virtual ~DistanceCalculator() {}
    virtual double calculate(const int i, const int j) = 0;
    // a clone method to make life easier when copying objects in each thread
    virtual DistanceCalculator* clone() const = 0;
    // helpers for distmat filler
    int xLimit() { return x_.length(); }
    int yLimit() { return y_.length(); }

protected:
    DistanceCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y)
        : dist_args_(DIST_ARGS)
        , x_(X)
        , y_(Y)
    { }

    int maxLength(const Rcpp::List& list, const bool is_multivariate) const {
        int max_len = 0;
        for (const SEXP& series : list) {
            if (is_multivariate) {
                Rcpp::NumericMatrix x(series);
                int this_len = x.nrow();
                if (this_len > max_len) max_len = this_len;
            }
            else {
                Rcpp::NumericVector x(series);
                int this_len = x.length();
                if (this_len > max_len) max_len = this_len;
            }
        }
        return max_len;
    }

    Rcpp::List dist_args_, x_, y_;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistanceCalculatorFactory
{
public:
    std::shared_ptr<DistanceCalculator> create(const SEXP& DIST, const SEXP& DIST_ARGS,
                                               const SEXP& X, const SEXP& Y);
    std::shared_ptr<DistanceCalculator> create(const std::string& dist, const SEXP& DIST_ARGS,
                                               const SEXP& X, const SEXP& Y);
};

// -------------------------------------------------------------------------------------------------
/* dtw_basic calculator */
// -------------------------------------------------------------------------------------------------
class DtwBasicCalculator : public DistanceCalculator
{
public:
    DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    DtwBasicCalculator* clone() const override
        { return new DtwBasicCalculator(*this); }

private:
    double calculate(const SEXP& X, const SEXP& Y);
    SEXP window_, norm_, step_, backtrack_, gcm_;
    bool is_multivariate_, normalize_;
};

// -------------------------------------------------------------------------------------------------
/* dtw_basic parallel calculator */
// -------------------------------------------------------------------------------------------------
class DtwBasicParallelCalculator : public DistanceCalculator
{
public:
    DtwBasicParallelCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    ~DtwBasicParallelCalculator();
    double calculate(const int i, const int j) override;
    DtwBasicParallelCalculator* clone() const override;

private:
    // method calculate
    double calculate(const RcppParallel::RVector<double>& x,
                     const RcppParallel::RVector<double>& y);
    double calculate(const RcppParallel::RMatrix<double>& x,
                     const RcppParallel::RMatrix<double>& y);
    // input parameters
    int window_;
    double norm_, step_;
    bool is_multivariate_;
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
/* lb_keogh calculator */
// -------------------------------------------------------------------------------------------------
class LbkCalculator : public DistanceCalculator
{
public:
    LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    LbkCalculator* clone() const override { return new LbkCalculator(*this); }

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
class LbiCalculator : public DistanceCalculator
{
public:
    LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    LbiCalculator* clone() const override { return new LbiCalculator(*this); }

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
class SdtwCalculator : public DistanceCalculator
{
public:
    SdtwCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    SdtwCalculator* clone() const override
        { return new SdtwCalculator(*this); }

private:
    double calculate(const SEXP& X, const SEXP& Y);
    SEXP gamma_, costmat_, mv_;
};

// -------------------------------------------------------------------------------------------------
/* gak calculator */
// -------------------------------------------------------------------------------------------------
class GakCalculator : public DistanceCalculator
{
public:
    GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    GakCalculator* clone() const override { return new GakCalculator(*this); }

private:
    double calculate(const SEXP& X, const SEXP& Y);
    SEXP sigma_, window_, logs_;
    bool is_multivariate_;
};

// -------------------------------------------------------------------------------------------------
/* sbd calculator */
// -------------------------------------------------------------------------------------------------
class SbdCalculator : public DistanceCalculator
{
public:
    SbdCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const int i, const int j) override;
    SbdCalculator* clone() const override { return new SbdCalculator(*this); }

private:
    double calculate(const arma::vec& x, const arma::vec& y,
                     const arma::cx_vec& fftx, const arma::cx_vec& ffty);
    arma::vec cc_seq_truncated_;
    Rcpp::List fftx_, ffty_;
    int fftlen_;
};

} // namespace dtwclust

#endif // DTWCLUST_DISTANCE_CALCULATORS_HPP_
