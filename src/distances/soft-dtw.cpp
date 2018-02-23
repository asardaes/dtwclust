#include "distances-details.h"

#include <algorithm> // std::max
#include <math.h> // exp, log, pow

#include <R.h> // R_PosInf

#include "../utils/SurrogateMatrix.h"

namespace dtwclust {

// =================================================================================================
/* point-wise squared Euclidean norm */
// =================================================================================================

double squared_euclidean(const double * const x, const double * const y,
                         const int i, const int j,
                         const int x_nrows, const int y_nrows, const int ncols)
{
    double d = 0;
    for (int k = 0; k < ncols; k++)
        d += pow(x[i + k * x_nrows] - y[j + k * y_nrows], 2);
    return d;
}

// =================================================================================================
/* soft min operator */
// =================================================================================================

double soft_min(double a, double b, double c, const double gamma)
{
    a /= -gamma;
    b /= -gamma;
    c /= -gamma;
    double max_val = std::max(std::max(a, b), c);
    double temp = 0;
    temp += exp(a - max_val);
    temp += exp(b - max_val);
    temp += exp(c - max_val);
    return -gamma * (log(temp) + max_val);
}

// =================================================================================================
/* thread-safe versions for the distance calculator */
// =================================================================================================

double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, SurrogateMatrix<double>& costmat)
{
    // initialize costmat values
    costmat[0] = 0;
    for (int i = 1; i < nx+2; i++) costmat(i,0) = R_PosInf;
    for (int j = 1; j < ny+2; j++) costmat(0,j) = R_PosInf;
    // compute distance
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            double point_norm = squared_euclidean(x, y, i-1, j-1, nx, ny, num_vars);
            costmat(i,j) = point_norm + soft_min(costmat(i-1,j),
                                                 costmat(i-1,j-1),
                                                 costmat(i,j-1),
                                                 gamma);
        }
    }
    return costmat(nx,ny);
}

double sdtw(const double * const x, const double * const y,
            const int nx, const int ny, const int num_vars,
            const double gamma, SurrogateMatrix<double>& costmat,
            SurrogateMatrix<double>& distmat)
{
    // initialize costmat values
    costmat[0] = 0;
    for (int i = 1; i < nx+2; i++) costmat(i,0) = R_PosInf;
    for (int j = 1; j < ny+2; j++) costmat(0,j) = R_PosInf;
    // compute distance
    for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
            double point_norm = squared_euclidean(x, y, i-1, j-1, nx, ny, num_vars);
            costmat(i,j) = point_norm + soft_min(costmat(i-1,j),
                                                 costmat(i-1,j-1),
                                                 costmat(i,j-1),
                                                 gamma);
            distmat(i-1,j-1) = point_norm;
        }
    }
    return costmat(nx,ny);
}

} // namespace dtwclust
