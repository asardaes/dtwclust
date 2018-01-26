#include "utils.h"

#include <algorithm> // std::fill
#include <vector>

namespace dtwclust {

KahanSummer::KahanSummer(double * const x, const int nrows, const int ncols)
    : x_(x)
    , nrows_(nrows)
    , ncols_(ncols)
    , c_(std::vector<double>(nrows * ncols))
    , y_(std::vector<double>(nrows * ncols))
    , t_(std::vector<double>(nrows * ncols))
{ }

void KahanSummer::reset() {
    std::fill(c_.begin(), c_.end(), 0);
}

void KahanSummer::add(const double value, const int i, const int j) {
    y_[d2s(i,j,nrows_)] = value - c_[d2s(i,j,nrows_)];
    t_[d2s(i,j,nrows_)] = x_[d2s(i,j,nrows_)] + y_[d2s(i,j,nrows_)];
    c_[d2s(i,j,nrows_)] = (t_[d2s(i,j,nrows_)] - x_[d2s(i,j,nrows_)]) - y_[d2s(i,j,nrows_)];
    x_[d2s(i,j,nrows_)] = t_[d2s(i,j,nrows_)];
}

} // namespace dtwclust
