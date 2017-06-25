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

RcppExport SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP DOTS);

RcppExport SEXP envelope(SEXP series, SEXP window);

RcppExport SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U);

RcppExport SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U);

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

// defined in envelope.cpp
void envelope_cpp(const Rcpp::NumericVector& array, const unsigned int width,
                  Rcpp::NumericVector& minvalues, Rcpp::NumericVector& maxvalues);

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
