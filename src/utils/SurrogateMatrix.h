#ifndef DTWCLUST_SURROGATEMATRIX_HPP_
#define DTWCLUST_SURROGATEMATRIX_HPP_

#include <algorithm> // std::fill

namespace dtwclust {

// a class to wrap raw pointers that should be treated as matrices
template<typename T>
class SurrogateMatrix
{
public:
    SurrogateMatrix() : x_(nullptr), own_x_(false) {}
    SurrogateMatrix(const int nrows, const int ncols, T * const x = nullptr)
        : x_(x ? x : new T[nrows * ncols])
        , nrows_(nrows)
        , ncols_(ncols)
        , own_x_(!x)
    { }

    SurrogateMatrix& operator=(const SurrogateMatrix&) = delete;

    ~SurrogateMatrix() {
        if (own_x_ && x_) delete[] x_;
    }

    SurrogateMatrix(const SurrogateMatrix& other) {
        nrows_ = other.nrows_;
        ncols_ = other.ncols_;
        own_x_ = other.own_x_;
        if (own_x_) {
            if (other.x_) {
                x_ = new T[nrows_ * ncols_];
                for (int id = 0; id < (nrows_ * ncols_); id++)
                    x_[id] = other.x_[id];
            }
            else {
                x_ = nullptr;
            }
        }
        else {
            x_ = other.x_;
        }
    }

    SurrogateMatrix(SurrogateMatrix&& other) {
        nrows_ = other.nrows_;
        ncols_ = other.ncols_;
        own_x_ = other.own_x_;
        x_ = other.x_;
        if (own_x_) other.x_ = nullptr;
    }

    SurrogateMatrix& operator=(SurrogateMatrix&& other) {
        nrows_ = other.nrows_;
        ncols_ = other.ncols_;
        own_x_ = other.own_x_;
        x_ = other.x_;
        if (own_x_) other.x_ = nullptr;
        return *this;
    }

    T& operator()(const int row, const int col) {
        return x_[row + col * nrows_];
    }

    const T operator()(const int row, const int col) const {
        return x_[row + col * nrows_];
    }

    T& operator[](const int id) {
        return x_[id];
    }

    const T operator[](const int id) const {
        return x_[id];
    }

    operator bool() const {
        return x_ ? true : false;
    }

    void fill(const T value) {
        std::fill(x_, x_ + (nrows_ * ncols_), value);
    }

private:
    T *x_;
    int nrows_, ncols_;
    bool own_x_;
};

} // namespace dtwclust

#endif // DTWCLUST_SURROGATEMATRIX_HPP_
