#include "details.h"

#include <cstddef> // std::size_t
#include <stdlib.h>
#include <math.h>

#include <R.h> // R_PosInf

#include "../utils/SurrogateMatrix.h"

namespace dtwclust {

// for cost matrix, in case of window constraint
#define NOT_VISITED -1.0

// for step matrix
#define UP 2.0
#define LEFT 1.0
#define DIAG 0.0

// double to single index, matrices are always vectors in R
std::size_t inline d2s(const std::size_t i, const std::size_t j, const std::size_t nx, const bool backtrack)
    __attribute__((always_inline));
std::size_t inline d2s(const std::size_t i, const std::size_t j, const std::size_t nx, const bool backtrack) {
    return backtrack ? (i + j * (nx + 1)) : ((i % 2) + j * 2);
}

// vector norm
double lnorm(const SurrogateMatrix<const double>& x, const SurrogateMatrix<const double>& y,
             const double norm, const std::size_t i, const std::size_t j)
{
    double res = 0;
    double temp;
    for (std::size_t k = 0; k < x.ncol(); k++) {
        temp = x(i,k) - y(j,k);

        if (norm == 1)
            temp = fabs(temp);
        else
            temp *= temp;

        res += temp;
    }

    return (norm == 1) ? res : sqrt(res);
}

// which direction to take in the cost matrix
int which_min(const double diag, const double left, const double up,
              const double step, const dtwclust_tuple_t local_cost,
              dtwclust_tuple_t * const tuple)
{
    tuple[0] = (diag == NOT_VISITED) ? R_PosInf : diag + step * local_cost;
    tuple[1] = (left == NOT_VISITED) ? R_PosInf : left + local_cost;
    tuple[2] = (up == NOT_VISITED) ? R_PosInf : up + local_cost;

    int direction = (tuple[1] < tuple[0]) ? 1 : 0;
    direction = (tuple[2] < tuple[direction]) ? 2 : direction;
    return direction;
}

// backtrack step matrix
int backtrack_steps(const SurrogateMatrix<double>& lcm,
                    SurrogateMatrix<int>& index1,
                    SurrogateMatrix<int>& index2,
                    const std::size_t nx,
                    const std::size_t ny)
{
    std::size_t i = nx - 1;
    std::size_t j = ny - 1;
    int path = 1;

    // always start at end of series
    index1[0] = nx;
    index2[0] = ny;
    while(!(i == 0 && j == 0)) {
        if (lcm[d2s(i, j, nx, true)] == 0) {
            i--;
            j--;
        }
        else if (lcm[d2s(i, j, nx, true)] == 1) {
            j--;
        }
        else if (lcm[d2s(i, j, nx, true)] == 2) {
            i--;
        }

        index1[path] = i + 1;
        index2[path] = j + 1;
        path++;
    }

    return path;
}

// the C code
double dtw_basic_c(SurrogateMatrix<double>& lcm,
                   const SurrogateMatrix<const double>& x,
                   const SurrogateMatrix<const double>& y,
                   const int w,
                   const double norm,
                   const double step,
                   const bool backtrack)
{
    std::size_t nx = x.nrow();
    std::size_t ny = y.nrow();
    dtwclust_tuple_t tuple[3];
    dtwclust_tuple_t local_cost;
    std::size_t i, j;
    int direction;

    // initialization (first row and first column)
    for (j = 0; j <= ny; j++) lcm[d2s(0, j, nx, backtrack)] = NOT_VISITED;
    for (i = 0; i <= (backtrack ? nx : 1); i++) lcm[d2s(i, 0, nx, backtrack)] = NOT_VISITED;

    // first value, must set here to avoid multiplying by step
    lcm[d2s(1, 1, nx, backtrack)] = lnorm(x, y, norm, 0, 0);
    if (norm == 2) lcm[d2s(1, 1, nx, backtrack)] *= lcm[d2s(1, 1, nx, backtrack)];

    // dynamic programming
    for (i = 1; i <= nx; i++) {
        std::size_t j1, j2;

        // adjust limits depending on window
        if (w == -1) {
            j1 = 1;
            j2 = ny;
        }
        else {
            // just use a double as temporary placeholder in case there are negatives
            double temp = ceil(static_cast<double>(i) * ny / nx - w);
            j1 = temp < 0 ? 0 : static_cast<std::size_t>(temp);
            j2 = static_cast<std::size_t>(floor(static_cast<double>(i) * ny / nx + w));
            j1 = j1 > 1 ? j1 : 1;
            j2 = j2 < ny ? j2 : ny;
        }

        for (j = 1; j <= ny; j++) {
            // very first value already set above
            if (i == 1 && j == 1) continue;

            // check if cell outside of window
            if (j < j1 || j > j2) {
                lcm[d2s(i, j, nx, backtrack)] = NOT_VISITED;
                continue;
            }

            local_cost = lnorm(x, y, norm, i-1, j-1);
            if (norm == 2) local_cost *= local_cost;

            direction = which_min(lcm[d2s(i-1, j-1, nx, backtrack)],
                                  lcm[d2s(i, j-1, nx, backtrack)],
                                  lcm[d2s(i-1, j, nx, backtrack)],
                                  step, local_cost, tuple);
            /*
             * I can use the same matrix to save both cost values and steps taken by shifting
             * the indices left and up for direction. Since the loop advances row-wise, the
             * appropriate values for the cost will be available, and the unnecessary ones are
             * replaced by steps along the way.
             */
            lcm[d2s(i, j, nx, backtrack)] = tuple[direction];
            if (backtrack) lcm[d2s(i-1, j-1, nx, backtrack)] = static_cast<double>(direction);
        }
    }
    return (norm == 1) ? lcm[d2s(nx, ny, nx, backtrack)] : sqrt(lcm[d2s(nx, ny, nx, backtrack)]);
}

// versions compatible with RcppParallel

double dtw_basic(SurrogateMatrix<double>& lcm,
                 const SurrogateMatrix<const double>& x,
                 const SurrogateMatrix<const double>& y,
                 const int window,
                 const double norm,
                 const double step,
                 const bool normalize,
                 const bool sqrt_dist)
{
    double d = dtw_basic_c(lcm, x, y, window, norm, step, false);
    if (!sqrt_dist) d *= d;
    if (normalize) d /= x.nrow() + y.nrow();
    return d;
}

double dtw_basic(SurrogateMatrix<double>& lcm,
                 const SurrogateMatrix<const double>& x,
                 const SurrogateMatrix<const double>& y,
                 const int window,
                 const double norm,
                 const double step,
                 const bool normalize,
                 const bool sqrt_dist,
                 SurrogateMatrix<int>& index1,
                 SurrogateMatrix<int>& index2,
                 int& path)
{
    double d = dtw_basic_c(lcm, x, y, window, norm, step, true);
    if (!sqrt_dist) d *= d;
    if (normalize) d /= x.nrow() + y.nrow();
    path = backtrack_steps(lcm, index1, index2, x.nrow(), y.nrow());
    return d;
}

} // namespace dtwclust
