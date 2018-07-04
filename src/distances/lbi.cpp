#include "details.h"

#include <cmath> // std::sqrt

#include "../utils/SurrogateMatrix.h"
#include "../utils/utils.h" // envelope_cpp, kahan_sum, id_t

namespace dtwclust {

// thread-safe
double lbi_core(const SurrogateMatrix<const double>& x,
                const SurrogateMatrix<const double>& y,
                const unsigned int window_size,
                const int p,
                const SurrogateMatrix<const double>& lower_envelope,
                const SurrogateMatrix<const double>& upper_envelope,
                SurrogateMatrix<double>& L2,
                SurrogateMatrix<double>& U2,
                SurrogateMatrix<double>& H,
                SurrogateMatrix<double>& LB)
{
    id_t length = x.nrow();
    double lb = 0;

    for (id_t i = 0; i < length; i++) {
        if (x[i] > upper_envelope[i]) {
            H[i] = upper_envelope[i];
            LB[i] = x[i] - upper_envelope[i];
        }
        else if (x[i] < lower_envelope[i]) {
            H[i] = lower_envelope[i];
            LB[i] = lower_envelope[i] - x[i];
        }
        else {
            H[i] = x[i];
            LB[i] = 0;
        }
        if (p > 1) LB[i] *= LB[i];
    }

    envelope_cpp(H, window_size * 2 + 1, L2, U2);
    double temp = 0;
    for (id_t i = 0; i < length; i++) {
        if (y[i] > U2[i])
            temp = y[i] - U2[i];
        else if (y[i] < L2[i])
            temp = L2[i] - y[i];
        else
            temp = 0;
        if (p > 1) temp *= temp;
        LB[i] += temp;
    }

    lb = kahan_sum(LB);
    if (p > 1) lb = std::sqrt(lb);
    return lb;
}

} // namespace dtwclust
