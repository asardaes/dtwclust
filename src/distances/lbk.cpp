#include "details.h"

#include <cmath> // std::sqrt

#include "../utils/SurrogateMatrix.h"
#include "../utils/utils.h" // kahan_sum, id_t

namespace dtwclust {

// thread-safe
double lbk_core(const SurrogateMatrix<const double>& x,
                const int p,
                const SurrogateMatrix<const double>& lower_envelope,
                const SurrogateMatrix<const double>& upper_envelope,
                SurrogateMatrix<double>& H)
{
    id_t length = x.nrow();
    double lb = 0;
    for (id_t i = 0; i < length; i++) {
        if (x[i] > upper_envelope[i])
            H[i] = x[i] - upper_envelope[i];
        else if (x[i] < lower_envelope[i])
            H[i] = lower_envelope[i] - x[i];
        else
            H[i] = 0;

        if (p > 1) H[i] *= H[i];
    }

    lb = kahan_sum(H);
    if (p > 1) lb = std::sqrt(lb);
    return lb;
}

} // namespace dtwclust
