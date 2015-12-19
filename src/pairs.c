#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rdefines.h>

/* ==================================================================================================================================================
 * All possible combinations in pairs
 * ==================================================================================================================================================
 */

void pairs_c(const int n, const int nrow, int *out, const int byrow)
{
     int i, j;
     int p = 0;

     if(byrow)
     {
          for(i = 1; i < n; i++)
          {
               for(j = i+1; j <= n; j++)
               {
                    out[p] = i;
                    out[p+nrow] = j;
                    p++;
               }
          }
     }
     else
     {
          for(j = 2; j <= n; j++)
          {
               for(i = 1; i < j; i++)
               {
                    out[p] = i;
                    out[p+nrow] = j;
                    p++;
               }
          }
     }
}

// the gateway function
SEXP pairs(SEXP L, SEXP byrow)
{
     int n = asInteger(L);
     int nrow = n * (n+1) / 2 - n;

     // allocate output integer vector
     SEXP ret = PROTECT(allocMatrix(INTSXP, nrow, 2));

     // get pointers to output objects
     int *out = INTEGER(ret);

     // dispatch to C function
     pairs_c(n, nrow, out, asLogical(byrow));

     // release protection
     UNPROTECT(1);

     // finish
     return ret;
}
