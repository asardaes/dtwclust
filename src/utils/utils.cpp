#include "utils.h"

#include <R.h>

namespace dtwclust {

// =================================================================================================
/* for Rcpp::Rcout */
// =================================================================================================

void Rflush()
{
    R_FlushConsole();
    R_ProcessEvents();
}

// =================================================================================================
/* helper kahan_sum */
// =================================================================================================

double kahan_sum(const double * const x, const int length)
{
    double sum = 0, c = 0;
    for (int i = 0; i < length; i++) {
        double y = x[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

// =================================================================================================
/* single to double indices for symmetric matrices without diagonal */
// =================================================================================================

void s2d(const int id, const int nrow, int& i, int& j)
{
    // check if it's the first column
    if (id < (nrow - 1)) {
        i = id + 1;
        j = 0;
        return;
    }
    // otherwise start at second column
    i = 2;
    j = 1;
    int start_id = nrow - 1;
    int end_id = nrow * 2 - 4;
    // j is ready after this while loop finishes
    while (!(id >= start_id && id <= end_id)) {
        start_id = end_id + 1;
        end_id = start_id + nrow - j - 3;
        i++;
        j++;
    }
    // while loop for i
    while (start_id < id) {
        i++;
        start_id++;
    }
}

} // namespace dtwclust
