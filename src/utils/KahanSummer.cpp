#include "utils.h"

#include <algorithm> // std::fill
#include <vector>

namespace dtwclust {

KahanSummer::KahanSummer(double * const x, const int nrows, const int ncols)
    : x_(x)
    , nrows_(nrows)
    , c_(std::vector<double>(nrows * ncols))
    , y_(std::vector<double>(nrows * ncols))
    , t_(std::vector<double>(nrows * ncols))
{ }

void KahanSummer::reset() {
    std::fill(c_.begin(), c_.end(), 0);
}

void KahanSummer::add(const double value, const int i, const int j) {
    int id = i + j * nrows_;
    y_[id] = value - c_[id];
    t_[id] = x_[id] + y_[id];
    c_[id] = (t_[id] - x_[id]) - y_[id];
    x_[id] = t_[id];
}

} // namespace dtwclust
