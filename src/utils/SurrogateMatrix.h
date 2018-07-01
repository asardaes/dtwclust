#ifndef DTWCLUST_SURROGATEMATRIX_HPP_
#define DTWCLUST_SURROGATEMATRIX_HPP_

#include <algorithm> // std::fill
#include <cstddef> // std::size_t

namespace dtwclust {

// a class to wrap raw pointers that should be treated as matrices
template<typename T>
class SurrogateMatrix
{
public:
    SurrogateMatrix() : x_(nullptr), own_x_(false) {}
    SurrogateMatrix(const std::size_t nrows, const std::size_t ncols)
        : x_(new T[nrows * ncols])
        , nrows_(nrows)
        , ncols_(ncols)
        , own_x_(true)
    { }
    SurrogateMatrix(const std::size_t nrows, const std::size_t ncols, T * const x)
        : x_(x)
        , nrows_(nrows)
        , ncols_(ncols)
        , own_x_(false)
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
                for (std::size_t id = 0; id < (nrows_ * ncols_); id++)
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

    std::size_t ncol() const {
        return ncols_;
    }

    std::size_t nrow() const {
        return nrows_;
    }

    T& operator()(const std::size_t row, const std::size_t col) {
        return x_[row + col * nrows_];
    }

    const T operator()(const std::size_t row, const std::size_t col) const {
        return x_[row + col * nrows_];
    }

    T& operator[](const std::size_t id) {
        return x_[id];
    }

    const T operator[](const std::size_t id) const {
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
    std::size_t nrows_, ncols_;
    bool own_x_;
};

} // namespace dtwclust

#endif // DTWCLUST_SURROGATEMATRIX_HPP_
