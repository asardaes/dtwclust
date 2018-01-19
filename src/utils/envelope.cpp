#include "utils.h"

#include <deque>

#include <RcppArmadillo.h>

namespace dtwclust {

/*
 * Lemire's streaming algorithm to compute warping envelope using no more than 3n comparisons
 *
 * Adapted from the code available at https://github.com/lemire/lbimproved/blob/master/dtw.h
 */

// thread-safe overload
void envelope_cpp(const double * const array, const int length, const unsigned int width,
                  double * const minvalues, double * const maxvalues)
{
    unsigned int constraint = (width - 1) / 2;
    unsigned int array_size = static_cast<unsigned int>(length);
    std::deque<int> maxfifo, minfifo;
    maxfifo.push_back(0);
    minfifo.push_back(0);
    for(unsigned int i = 1; i < array_size; ++i) {
        if(i >= constraint + 1) {
            maxvalues[i - constraint - 1] = array[maxfifo.front()];
            minvalues[i - constraint - 1] = array[minfifo.front()];
        }
        if(array[i] > array[i - 1]) { //overshoot
            maxfifo.pop_back();
            while(maxfifo.size() > 0) {
                if(array[i] <= array[maxfifo.back()]) break;
                maxfifo.pop_back();
            }
        }
        else {
            minfifo.pop_back();
            while(minfifo.size() > 0) {
                if(array[i] >= array[minfifo.back()]) break;
                minfifo.pop_back();
            }
        }
        maxfifo.push_back(i);
        minfifo.push_back(i);
        if(i == width + maxfifo.front())
            maxfifo.pop_front();
        else if(i == width + minfifo.front())
            minfifo.pop_front();
    }
    for(unsigned int i = array_size; i <= array_size + constraint; ++i) {
        maxvalues[i - constraint - 1] = array[maxfifo.front()];
        minvalues[i - constraint - 1] = array[minfifo.front()];
        if(i - maxfifo.front() >= width) maxfifo.pop_front();
        if(i - minfifo.front() >= width) minfifo.pop_front();
    }
}

// non-thread-safe
void envelope_cpp(const Rcpp::NumericVector& array, const unsigned int width,
                  Rcpp::NumericVector& minvalues, Rcpp::NumericVector& maxvalues)
{
    int length = array.length();
    envelope_cpp(&array[0], length, width, &minvalues[0], &maxvalues[0]);
}

// gateway
RcppExport SEXP envelope(SEXP series, SEXP window) {
    BEGIN_RCPP
    Rcpp::NumericVector x(series);
    Rcpp::NumericVector L(x.size()), U(x.size());
    envelope_cpp(x, Rcpp::as<unsigned int>(window), L, U);
    Rcpp::List ret;
    ret["lower"] = L;
    ret["upper"] = U;
    return(ret);
    END_RCPP
}

} // namespace dtwclust
