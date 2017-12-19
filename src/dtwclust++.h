#include <Rcpp.h>
#include <algorithm> // stable_sort
#include <numeric> // iota
#include <utility> // swap
#include <vector>

#ifndef _DTWCLUST_HPP
#define _DTWCLUST_HPP

namespace dtwclust {

// =================================================================================================
/* Functions that can be called from R */
// =================================================================================================

RcppExport SEXP SparseDistmatIndices__new(SEXP num_rows);
RcppExport SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric);

RcppExport SEXP dba(SEXP X, SEXP centroid,
                    SEXP max_iter, SEXP delta, SEXP trace,
                    SEXP multivariate, SEXP mv_ver, SEXP DOTS);

RcppExport SEXP dtwb_loop(SEXP D, SEXP X, SEXP Y,
                          SEXP SYMMETRIC, SEXP PAIRWISE,
                          SEXP BIGMAT, SEXP NORMALIZE, SEXP MULTIVARIATE,
                          SEXP DISTARGS, SEXP ENDPOINTS);

RcppExport SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS);

RcppExport SEXP envelope(SEXP series, SEXP window);

RcppExport SEXP force_lb_symmetry(SEXP X);

RcppExport SEXP gak_loop(SEXP D, SEXP X, SEXP Y,
                         SEXP SYMMETRIC, SEXP PAIRWISE,
                         SEXP BIGMAT, SEXP MULTIVARIATE,
                         SEXP DISTARGS, SEXP ENDPOINTS);

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbk_loop(SEXP D, SEXP X, SEXP L, SEXP U,
                         SEXP PAIRWISE, SEXP BIGMAT,
                         SEXP P, SEXP LEN, SEXP ENDPOINTS);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbi_loop(SEXP D, SEXP X, SEXP Y, SEXP L, SEXP U,
                         SEXP PAIRWISE, SEXP BIGMAT,
                         SEXP P, SEXP WINDOW, SEXP LEN, SEXP ENDPOINTS);

RcppExport SEXP sbd_loop(SEXP D, SEXP X, SEXP Y, SEXP FFTX, SEXP FFTY,
                         SEXP FFTLEN, SEXP SYMMETRIC, SEXP PAIRWISE, SEXP ENDPOINTS, SEXP BIGMAT);

RcppExport SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID, SEXP GAMMA, SEXP WEIGHTS, SEXP MV,
                          SEXP COSTMAT, SEXP DISTMAT, SEXP EM);

RcppExport SEXP sdtw_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                          SEXP SYMMETRIC, SEXP PAIRWISE,
                          SEXP BIGMAT, SEXP MV,
                          SEXP ENDPOINTS);

RcppExport SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP DISTMAT, SEXP MV);

RcppExport SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST);

// =================================================================================================
/* Internal functions that are called in more than one C++ file */
// =================================================================================================

// defined in utils.cpp
void Rflush();

// defined in utils.cpp
double dtwb(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& dots);
double dtwb(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y, const Rcpp::List& dots);

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

// defined in utils.cpp
double gak(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& dots);
double gak(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y, const Rcpp::List& dots);

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
        // while vOrder[i] is not yet in place
        // every swap places at least one element in its proper place
        while (order[i] != order[order[i]]) {
            std::swap(v[order[i]], v[order[order[i]]]);
            std::swap(order[i], order[order[i]]);
        }
    }
}

} // namespace dtwclust

#endif
