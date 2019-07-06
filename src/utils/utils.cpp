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
void Rflush() {
    R_FlushConsole();
    R_ProcessEvents();
}

/* helper kahan_sum */
double kahan_sum(const SurrogateMatrix<double>& x) {
    double sum = 0, c = 0;
    for (id_t i = 0; i < x.nrow(); i++) {
        double y = x[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

/**
 * single to double indices for symmetric matrices without diagonal
 * https://stackoverflow.com/a/27088560/5793905
 * https://math.stackexchange.com/questions/646117/how-to-find-a-function-mapping-matrix-indices
 */
void s2d(const id_t id, const id_t nrow, id_t& i, id_t& j) {
    // tmp to avoid ambiguity
    double tmp = static_cast<double>(-8 * id + 4 * nrow * (nrow - 1) - 7);
    j = nrow - 2 - static_cast<id_t>(sqrt(tmp) / 2 - 0.5);
    i = id + j + 1 - nrow * (nrow - 1) / 2 + (nrow - j) * ((nrow - j) - 1) / 2;
}

} // namespace dtwclust
