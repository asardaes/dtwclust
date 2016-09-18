#include <stdlib.h>
#include <math.h>
#include "shared.h"

double lnorm(const double *x, const double *y, const double norm,
             const int nx, const int ny, const int dim,
             const int i, const int j)
{
     double res = 0;
     for (int k = 0; k < dim; k++) res += pow(fabs(x[i + nx*k] - y[j + ny*k]), norm);
     res = pow(res, 1/norm);
     return res;
}
