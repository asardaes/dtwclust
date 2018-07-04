#include "utils.h"

#include <deque>

#include "SurrogateMatrix.h"

namespace dtwclust {

/*
 * Lemire's streaming algorithm to compute warping envelope using no more than 3n comparisons
 *
 * Adapted from the code available at https://github.com/lemire/lbimproved/blob/master/dtw.h
 */

// thread-safe
void envelope_cpp(const SurrogateMatrix<double>& array, const unsigned int width,
                  SurrogateMatrix<double>& minvalues, SurrogateMatrix<double>& maxvalues)
{
    unsigned int constraint = (width - 1) / 2;
    id_t array_size = array.nrow();
    std::deque<int> maxfifo, minfifo;
    maxfifo.push_back(0);
    minfifo.push_back(0);
    for(id_t i = 1; i < array_size; ++i) {
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
    for(id_t i = array_size; i <= array_size + constraint; ++i) {
        maxvalues[i - constraint - 1] = array[maxfifo.front()];
        minvalues[i - constraint - 1] = array[minfifo.front()];
        if(i - maxfifo.front() >= width) maxfifo.pop_front();
        if(i - minfifo.front() >= width) minfifo.pop_front();
    }
}

} // namespace dtwclust
