/*
 * Lemire's streaming algorithm to compute warping envelop using no more than 3n comparisons
 *
 * Adapted from the code available at https://github.com/lemire/lbimproved/blob/master/dtw.h
 */

#include <Rcpp.h>
#include <deque>

using namespace std;
using namespace Rcpp;

void envelop_cpp(const NumericVector& array, unsigned int width,
                 NumericVector& minvalues, NumericVector& maxvalues) {

     unsigned int constraint = (width - 1) / 2;
     unsigned int array_size = static_cast<unsigned int>(array.size());
     deque<int> maxfifo, minfifo;

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

          } else {
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

     for(unsigned int i = array.size(); i <= array_size + constraint; ++i) {
          maxvalues[i - constraint - 1] = array[maxfifo.front()];
          minvalues[i - constraint - 1] = array[minfifo.front()];

          if(i - maxfifo.front() >= width) maxfifo.pop_front();
          if(i - minfifo.front() >= width) minfifo.pop_front();
     }
}

RcppExport SEXP envelop (SEXP series, SEXP window) {
     NumericVector x(series);
     NumericVector L(x.size());
     NumericVector U(x.size());

     envelop_cpp(x, as<unsigned int>(window), L, U);

     List ret;
     ret["min"] = L;
     ret["max"] = U;

     return(ret);
}
