#include "distances-details.h"

#include <stdlib.h>
#include <math.h>
#include <type_traits> // conditional

#include <R.h> // R_PosInf

namespace dtwclust {

// for cost matrix, in case of window constraint
#define NOT_VISITED -1.0

// for step matrix
#define UP 2.0
#define LEFT 1.0
#define DIAG 0.0

// double to single index, matrices are always vectors in R
int inline d2s(int const i, int const j, int const nx, int const backtrack)
    __attribute__((always_inline));
int inline d2s(int const i, int const j, int const nx, int const backtrack) {
    return backtrack ? (i + j * (nx + 1)) : ((i % 2) + j * 2);
}

// vector norm
double lnorm(double const * const x, double const * const y, double const norm,
             int const nx, int const ny, int const num_var,
             int const i, int const j)
{
    double res = 0;
    double temp;
    for (int k = 0; k < num_var; k++) {
        temp = x[i + nx * k] - y[j + ny * k];
        if (norm == 1)
            temp = fabs(temp);
        else
            temp *= temp;
        res += temp;
    }
    return (norm == 1) ? res : sqrt(res);
}

// which direction to take in the cost matrix
int which_min(double const diag, double const left, double const up,
              double const step, dtwclust_tuple_t const local_cost, dtwclust_tuple_t * const tuple)
{
    // DIAG, LEFT, UP
    tuple[0] = (diag == NOT_VISITED) ? R_PosInf : diag + step * local_cost;
    tuple[1] = (left == NOT_VISITED) ? R_PosInf : left + local_cost;
    tuple[2] = (up == NOT_VISITED) ? R_PosInf : up + local_cost;

    int direction = (tuple[1] < tuple[0]) ? 1 : 0;
    direction = (tuple[2] < tuple[direction]) ? 2 : direction;
    return direction;
}

// backtrack step matrix
int backtrack_steps(double const * const D,
                    int const nx, int const ny,
                    int *index1, int *index2)
{
    int i = nx - 1;
    int j = ny - 1;
    int path = 1;

    // always start at end of series
    index1[0] = nx;
    index2[0] = ny;
    while(!(i == 0 && j == 0)) {
        if (D[d2s(i, j, nx, 1)] == 0) {
            i--;
            j--;
        }
        else if (D[d2s(i, j, nx, 1)] == 1) {
            j--;
        }
        else if (D[d2s(i, j, nx, 1)] == 2) {
            i--;
        }
        index1[path] = i + 1;
        index2[path] = j + 1;
        path++;
    }
    return path;
}

// the C code
double dtw_basic_c(double * const D, dtwclust_tuple_t * const tuple,
                   double const * const x, double const * const y, int const w,
                   int const nx, int const ny, int const num_var,
                   double const norm, double const step,
                   int const backtrack)
{
    int i, j, direction;
    dtwclust_tuple_t local_cost;
    // initialization (first row and first column)
    for (j = 0; j <= ny; j++) D[d2s(0, j, nx, backtrack)] = NOT_VISITED;
    for (i = 0; i <= (backtrack ? nx : 1); i++) D[d2s(i, 0, nx, backtrack)] = NOT_VISITED;
    // first value, must set here to avoid multiplying by step
    D[d2s(1, 1, nx, backtrack)] = lnorm(x, y, norm, nx, ny, num_var, 0, 0);
    if (norm == 2) D[d2s(1, 1, nx, backtrack)] *= D[d2s(1, 1, nx, backtrack)];
    // dynamic programming
    for (i = 1; i <= nx; i++) {
        int j1, j2;
        // adjust limits depending on window
        if (w == -1) {
            j1 = 1;
            j2 = ny;
        }
        else {
            j1 = ceil((double) i * ny / nx - w);
            j2 = floor((double) i * ny / nx + w);
            j1 = j1 > 1 ? j1 : 1;
            j2 = j2 < ny ? j2 : ny;
        }
        for (j = 1; j <= ny; j++) {
            // very first value already set above
            if (i == 1 && j == 1) continue;
            // check if cell outside of window
            if (j < j1 || j > j2) {
                D[d2s(i, j, nx, backtrack)] = NOT_VISITED;
                continue;
            }
            local_cost = lnorm(x, y, norm, nx, ny, num_var, i-1, j-1);
            if (norm == 2) local_cost *= local_cost;
            // set the value of 'direction'
            direction = which_min(D[d2s(i-1, j-1, nx, backtrack)], D[d2s(i, j-1, nx, backtrack)],
                                  D[d2s(i-1, j, nx, backtrack)], step, local_cost, tuple);
            /*
            * I can use the same matrix to save both cost values and steps taken by shifting
            * the indices left and up for direction. Since the loop advances row-wise, the
            * appropriate values for the cost will be available, and the unnecessary ones are
            * replaced by steps along the way.
            */
            D[d2s(i, j, nx, backtrack)] = tuple[direction];
            if (backtrack) D[d2s(i-1, j-1, nx, backtrack)] = (double) direction;
        }
    }
    return (norm == 1) ? D[d2s(nx, ny, nx, backtrack)] : sqrt(D[d2s(nx, ny, nx, backtrack)]);
}

// a version compatible with RcppParallel
double dtw_basic_par(double const * const x, double const * const y,
                     int const nx, int const ny, int const num_var,
                     int const window, double const norm, double const step, int const normalize,
                     double * const distmat, int const backtrack,
                     int * const index1, int * const index2, int * const path)
{
    dtwclust_tuple_t tuple[3];
    double d = dtw_basic_c(distmat, tuple, x, y, window, nx, ny, num_var, norm, step, backtrack);
    if (normalize) d /= nx + ny;
    if (backtrack) *path = backtrack_steps(distmat, nx, ny, index1, index2);
    return d;
}

} // namespace dtwclust
