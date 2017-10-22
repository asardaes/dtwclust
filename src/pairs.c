#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "dtwclust.h"

/* =================================================================================================
 * All possible combinations in pairs
 * =================================================================================================
 */
void pairs_c(const int n, const int nrow, int *out)
{
    int i, j;
    int p = 0;
    for(j = 1; j < n; j++)
    {
        for(i = j+1; i <= n; i++)
        {
            out[p] = i;
            out[p+nrow] = j;
            p++;
        }
    }
}

// the gateway function
SEXP pairs(SEXP L)
{
    int n = asInteger(L);
    int nrow = n * (n+1) / 2 - n;

    // allocate output integer vector
    SEXP ret = PROTECT(allocMatrix(INTSXP, nrow, 2));

    // dispatch to C function
    pairs_c(n, nrow, INTEGER(ret));

    // release protection
    UNPROTECT(1);

    // finish
    return ret;
}
