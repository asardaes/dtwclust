#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

/* ==================================================================================================================================================
 * Update local cost matrix for DBA
 * ==================================================================================================================================================
 */

void update_lcm_c(double *lcm, const double *x, const double *c, const int sqr, const int I, const int J)
{
     int i, j;

     for(i = 0; i < I; i++)
     {
          for(j = 0; j < J; j++)
          {
               if(sqr)
                    lcm[i + I*j] = pow(x[i] - c[j], 2);
               else
                    lcm[i + I*j] = fabs(x[i] - c[j]);
          }
     }
}

// the gateway function
SEXP update_lcm(SEXP lcm, SEXP x, SEXP center, SEXP square)
{
     int I, J;

     I = LENGTH(x);
     J = LENGTH(center);

     if(LENGTH(lcm) != I*J)
          error("Invalid dimensions for local cost matrix when computing DBA");

     // dispatch to C function
     update_lcm_c(REAL(lcm), REAL(x), REAL(center), asLogical(square), I, J);

     // finish
     return R_NilValue;
}
