#include "details.h"

#include <algorithm> // std::max
#include <math.h> // exp, log, pow

#include <R.h> // R_PosInf

#include "../utils/SurrogateMatrix.h"
#include "../utils/utils.h" // id_t

namespace dtwclust {

// =================================================================================================
/* point-wise squared Euclidean norm */
// =================================================================================================

double squared_euclidean(const SurrogateMatrix<const double>& x,
                         const SurrogateMatrix<const double>& y,
                         const id_t i, const id_t j)
{
    double d = 0;
    for (id_t k = 0; k < x.ncol(); k++)
        d += pow(x(i,k) - y(j,k), 2);
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

double sdtw(const SurrogateMatrix<const double>& x, const SurrogateMatrix<const double>& y,
            const double gamma, SurrogateMatrix<double>& costmat)
{
    SurrogateMatrix<double> dummy_distmat;
    return sdtw(x, y, gamma, costmat, dummy_distmat);
}

double sdtw(const SurrogateMatrix<const double>& x, const SurrogateMatrix<const double>& y,
            const double gamma, SurrogateMatrix<double>& costmat, SurrogateMatrix<double>& distmat)
{
    id_t nx = x.nrow(), ny = y.nrow();
    // initialize costmat values
    costmat[0] = 0;
    for (id_t i = 1; i < nx+2; i++) costmat(i,0) = R_PosInf;
    for (id_t j = 1; j < ny+2; j++) costmat(0,j) = R_PosInf;
    // compute distance
    for (id_t i = 1; i <= nx; i++) {
        for (id_t j = 1; j <= ny; j++) {
            double point_norm = squared_euclidean(x, y, i-1, j-1);
            costmat(i,j) = point_norm + soft_min(costmat(i-1,j),
                                                 costmat(i-1,j-1),
                                                 costmat(i,j-1),
                                                 gamma);
            if (distmat) distmat(i-1,j-1) = point_norm;
        }
    }
    return costmat(nx,ny);
}

} // namespace dtwclust
