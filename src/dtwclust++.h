#include <RcppArmadillo.h> // sbd uses fft
#include <algorithm> // stable_sort
#include <memory> // *_ptr
#include <numeric> // iota
#include <string>
#include <utility> // swap
#include <vector>

#ifndef _DTWCLUST_HPP_
#define _DTWCLUST_HPP_

namespace dtwclust {

// =================================================================================================
/* Functions that can be called from R */
// =================================================================================================

RcppExport SEXP SparseDistmatIndices__new(SEXP num_rows);
RcppExport SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

RcppExport SEXP dba(SEXP X, SEXP centroid,
                    SEXP max_iter, SEXP delta, SEXP trace,
                    SEXP multivariate, SEXP mv_ver, SEXP DOTS);

RcppExport SEXP dtwb_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                          SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS);

RcppExport SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS);

RcppExport SEXP envelope(SEXP series, SEXP window);

RcppExport SEXP force_lb_symmetry(SEXP X);

RcppExport SEXP gak_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                         SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS);

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbk_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                         SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbi_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                         SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS);

RcppExport SEXP sbd_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                         SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS);

RcppExport SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID, SEXP GAMMA, SEXP WEIGHTS, SEXP MV,
                          SEXP COSTMAT, SEXP DISTMAT, SEXP EM);

RcppExport SEXP sdtw_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                          SEXP SYMMETRIC, SEXP PAIRWISE, SEXP BIGMAT, SEXP ENDPOINTS);

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP DISTMAT, SEXP MV);

RcppExport SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST);

// =================================================================================================
/* Enums */
// =================================================================================================

enum class Distance {
    DTW_BASIC,
    LBK,
    LBI,
    SDTW,
    GAK,
    SBD
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
    virtual ~DistanceCalculator() {};
    virtual double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
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
    std::shared_ptr<DistanceCalculator>
    createCalculator(enum Distance distance, const SEXP& DIST_ARGS);
};

// -------------------------------------------------------------------------------------------------
/* dtw_basic calculator */
// -------------------------------------------------------------------------------------------------
class DtwBasicDistanceCalculator : public DistanceCalculator
{
public:
    DtwBasicDistanceCalculator(const SEXP& DIST_ARGS);
    double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) override;
private:
    double calculateDistance(const SEXP& X, const SEXP& Y);
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
    double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) override;
private:
    double calculateDistance(const Rcpp::NumericVector& x,
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
    double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) override;
private:
    double calculateDistance(const Rcpp::NumericVector& x,
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
    double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) override;
private:
    double calculateDistance(const SEXP& X, const SEXP& Y);
    SEXP gamma_, costmat_, mv_;
};

// -------------------------------------------------------------------------------------------------
/* gak calculator */
// -------------------------------------------------------------------------------------------------
class GakDistanceCalculator : public DistanceCalculator
{
public:
    GakDistanceCalculator(const SEXP& DIST_ARGS);
    double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) override;
private:
    double calculateDistance(const SEXP& X, const SEXP& Y);
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
    double calculateDistance(const Rcpp::List& X, const Rcpp::List& Y,
                             const int i, const int j) override;
private:
    double calculateDistance(const arma::vec& x, const arma::vec& y,
                             const arma::cx_vec& fftx, const arma::cx_vec& ffty);
    arma::vec cc_seq_truncated_;
    Rcpp::List fftx_, ffty_;
    int fftlen_;
};

// =================================================================================================
/* DistmatFillStrategy (base + concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distmat fill strategy */
// -------------------------------------------------------------------------------------------------
class DistmatFillStrategy
{
public:
    virtual ~DistmatFillStrategy() {};
    virtual void fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                             const std::shared_ptr<DistanceCalculator>& dist_calculator,
                             const SEXP& ENDPOINTS, const bool is_bigmat) const = 0;
};

// -------------------------------------------------------------------------------------------------
/* pairwise distmat fill strategy */
// -------------------------------------------------------------------------------------------------
class PairwiseDistmatFill : public DistmatFillStrategy
{
public:
    void fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                     const std::shared_ptr<DistanceCalculator>& dist_calculator,
                     const SEXP& ENDPOINTS, const bool is_bigmat) const override;
};

// -------------------------------------------------------------------------------------------------
/* symmetric distmat fill strategy */
// -------------------------------------------------------------------------------------------------
class SymmetricDistmatFill : public DistmatFillStrategy
{
public:
    void fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                     const std::shared_ptr<DistanceCalculator>& dist_calculator,
                     const SEXP& ENDPOINTS, const bool is_bigmat) const override;
};

// -------------------------------------------------------------------------------------------------
/* general distmat fill strategy */
// -------------------------------------------------------------------------------------------------
class GeneralDistmatFill : public DistmatFillStrategy
{
public:
    void fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                     const std::shared_ptr<DistanceCalculator>& dist_calculator,
                     const SEXP& ENDPOINTS, const bool is_bigmat) const override;
};

// =================================================================================================
/* DistmatFiller */
// =================================================================================================

class DistmatFiller
{
public:
    DistmatFiller(const SEXP& IS_BIGMAT, const SEXP& ENDPOINTS,
                  enum Distance distance, const SEXP& DIST_ARGS);

    void chooseFillStrategy(const bool pairwise, const bool symmetric);
    void fillDistmat(const SEXP& D, const SEXP& X, const SEXP& Y) const;

private:
    std::unique_ptr<DistmatFillStrategy> fill_strategy_;
    std::shared_ptr<DistanceCalculator> dist_calculator_;
    SEXP ENDPOINTS_;
    bool is_bigmat_;
};

// =================================================================================================
/* Internal functions that are called in more than one C++ file */
// =================================================================================================

// defined in utils.cpp
void Rflush();

// defined in utils.cpp
void fill_distmat(DistmatFiller& distmat_filler,
                  const SEXP& D, const SEXP& X, const SEXP& Y,
                  const bool pairwise, const bool symmetric);

// defined in envelope.cpp
void envelope_cpp(const Rcpp::NumericVector& array, const unsigned int width,
                  Rcpp::NumericVector& minvalues, Rcpp::NumericVector& maxvalues);

// defined in lbs.cpp
double lbk_core(const Rcpp::NumericVector& x, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& H);

double lbi_core(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y,
                const unsigned int window_size, const int p,
                const Rcpp::NumericVector& lower_envelope,
                const Rcpp::NumericVector& upper_envelope,
                Rcpp::NumericVector& L2,
                Rcpp::NumericVector& U2,
                Rcpp::NumericVector& H,
                Rcpp::NumericVector& LB);

// =================================================================================================
/* Templates */
// =================================================================================================

// see https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> stable_sort_ind(const std::vector<T>& v, const bool decreasing)
{
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    if (decreasing)
        std::stable_sort(idx.begin(), idx.end(),
                         [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
    else
        std::stable_sort(idx.begin(), idx.end(),
                         [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}

// see https://stackoverflow.com/a/22183350/5793905
template <typename T>
void reorder(std::vector<T>& v, std::vector<size_t>& order)
{
    if (v.size() != order.size()) Rcpp::stop("Cannot reorder with vectors of different sizes");

    // for all elements to put in place
    for (size_t i = 0; i < v.size(); ++i) {
        // while order[i] is not yet in place
        // every swap places at least one element in its proper place
        while (order[i] != order[order[i]]) {
            std::swap(v[order[i]], v[order[order[i]]]);
            std::swap(order[i], order[order[i]]);
        }
    }
}

} // namespace dtwclust

#endif // _DTWCLUST_HPP_
