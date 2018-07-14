#include "utils.h"

#include <R.h>

#include "SurrogateMatrix.h"

namespace dtwclust {

/* grain parameter for multi-threading */
int get_grain(const int n, const int num_threads) {
    int grain1 = n / num_threads / 10;
    int grain2 = n / 100;
    int grain = (grain1 < grain2) ? grain1 : grain2;
    return (grain < DTWCLUST_MIN_GRAIN) ? DTWCLUST_MIN_GRAIN : grain;
}

/* for Rcpp::Rcout */
void Rflush()
{
    R_FlushConsole();
    R_ProcessEvents();
}

/* helper kahan_sum */
double kahan_sum(const SurrogateMatrix<double>& x)
{
    double sum = 0, c = 0;
    for (id_t i = 0; i < x.nrow(); i++) {
        double y = x[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

/* single to double indices for symmetric matrices without diagonal */
void s2d(const id_t id, const id_t nrow, id_t& i, id_t& j) {
    // check if it's the first column
    if (id < (nrow - 1)) {
        i = id + 1;
        j = 0;
        return;
    }
    // otherwise start at second column
    i = 2;
    j = 1;
    id_t start_id = nrow - 1;
    id_t end_id = nrow * 2 - 4;
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
