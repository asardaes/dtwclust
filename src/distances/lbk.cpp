#include "distances-details.h"

#include <cmath> // std::sqrt

#include "../utils/utils.h" // kahan_sum

namespace dtwclust {

// thread-safe
double lbk_core(const double * const x, const int length, const int p,
                const double * const lower_envelope, const double * const upper_envelope,
                double * const H)
{
    double lb = 0;
    for (int i = 0; i < length; i++) {
        if (x[i] > upper_envelope[i])
            H[i] = x[i] - upper_envelope[i];
        else if (x[i] < lower_envelope[i])
            H[i] = lower_envelope[i] - x[i];
        else
            H[i] = 0;
        if (p > 1) H[i] *= H[i];
    }
    lb = kahan_sum(H, length);
    if (p > 1) lb = std::sqrt(lb);
    return lb;
}

} // namespace dtwclust
