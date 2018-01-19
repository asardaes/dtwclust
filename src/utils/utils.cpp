#include "utils.h"

#include <RcppArmadillo.h>

namespace dtwclust {

// =================================================================================================
/* Force symmetry helper */
// =================================================================================================

RcppExport SEXP force_lb_symmetry(SEXP X)
{
    BEGIN_RCPP
    Rcpp::NumericMatrix matrix(X);
    for (int i = 1; i < matrix.nrow(); i++) {
        Rcpp::checkUserInterrupt();
        for (int j = 0; j < i; j++) {
            double lb1 = matrix(i,j);
            double lb2 = matrix(j,i);
            if (lb1 > lb2)
                matrix(j,i) = lb1;
            else
                matrix(i,j) = lb2;
        }
    }
    return R_NilValue;
    END_RCPP
}

// =================================================================================================
/* all possible combinations in pairs */
// =================================================================================================

void pairs_c(const int n, const int nrow, int *out)
{
    int i, j;
    int p = 0;
    for(j = 1; j < n; j++)
    {
        for(i = j + 1; i <= n; i++)
        {
            out[p] = i;
            out[p+nrow] = j;
            p++;
        }
    }
}

// the gateway function
RcppExport SEXP pairs(SEXP L)
{
    int n = Rf_asInteger(L);
    int nrow = n * (n+1) / 2 - n;

    // allocate output integer vector
    SEXP ret = PROTECT(Rf_allocMatrix(INTSXP, nrow, 2));

    // dispatch to C function
    pairs_c(n, nrow, INTEGER(ret));

    // release protection
    UNPROTECT(1);

    // finish
    return ret;
}

// =================================================================================================
/* assign existing names to existing vector */
// =================================================================================================

RcppExport SEXP setnames_inplace(SEXP vec, SEXP names) {
    Rf_setAttrib(vec, R_NamesSymbol, names);
    return R_NilValue;
}

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
