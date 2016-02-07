#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <float.h>

/* ==================================================================================================================================================
 * Running extremes
 *
 * Adapted from the runmin/runmax functions that can be found in the R package "caTools" by Jarek Tuszynski.
 * ==================================================================================================================================================
 */

void runminmax_c(double *in, double *outMin, double *outMax, const int n, const int m)
{
     /* full-blown version with NaN's and edge calculation */
     int i, j, k2;
     double ptOutMin, ptOutMax, Min, Max, CSTmin = DBL_MAX, CSTmax = -DBL_MAX;
     double *inMin = in, *inMax = in; // input vector doesn't change, so pointers can point to the same
     double NaN = R_NaN;

     k2  = m>>1;               /* right half of window size */

     /* --- step 1 - find min/max of elements 0:(k2-1) */
     Min = CSTmin;                  /* we need to calculate  initial 'Min' */
     Max = CSTmax;                /* we need to calculate  initial 'Max' */

     for(i=0; i<k2; i++){
          Min = fmin(Min,inMin[i]);  /* find minimum over a window of length k2 */
     Max = fmax(Max,inMax[i]);  /* find maximum over a window of length k2 */
     }

     /* --- step 2 - left edge - start expanding the moving window to the right */
     for(i=k2; i<m-1; i++) {
          Min=fmin(Min,inMin[i]);     /* cumulative min */
     *(outMin++) = (Min==CSTmin ? NaN : Min); /* save 'Min' and move window */

     Max=fmax(Max,inMax[i]);     /* cumulative max */
     *(outMax++) = (Max==CSTmax ? NaN : Max); /* save 'Max' and move window */
     }

     /* --- step 3 - the inner section - window of constant size is moving  */
     ptOutMin=CSTmin;
     ptOutMax=CSTmax;

     for(i=m-1; i<n; i++) {
          if(ptOutMin==Min) {        /* if point comining outMin of the window was window's min than ... */
     Min=CSTmin;              /* we need to recalculate 'Min' */
     for(j=0; j<m; j++)
          Min=fmin(Min,inMin[j]); /* find minimum over a window of length m */
          } else                  /* if point comining outMin of the window was NOT window min than min of ... */
     Min=fmin(Min,inMin[m-1]); /* ... window's first m-1 points is still 'Min', so we have to add a single point */
     ptOutMin = *(inMin++);        /* store point comming outMin of the window for future use and move window */
     *(outMin++) = (Min==CSTmin ? NaN : Min); /* save 'Min' and move window */

     if(ptOutMax==Max) {        /* if point comaxing outMax of the window was window's max than ... */
     Max=CSTmax;              /* we need to recalculate 'Max' */
     for(j=0; j<m; j++)
          Max=fmax(Max,inMax[j]); /* find maximum over a window of length m */
     } else                  /* if point comining outMax of the window was NOT window max than max of ... */
     Max=fmax(Max,inMax[m-1]); /* ... window's first m-1 points is still 'Max', so we have to add a single point */
     ptOutMax = *(inMax++);        /* store point comming outMax of the window for future use and move window */
     *(outMax++) = (Max==CSTmax ? NaN : Max); /* save 'Max' and move window */
     }

     /* --- step 4 - right edge - right side reached the end and left is shrinking  */
     for(i=0; i<k2; i++) {
          if(ptOutMin==Min) {        /* if point comining outMin of the window was window's extreme than ... */
     Min=CSTmin;              /* we need to recalculate 'Min' */
     for(j=0; j<m-i-1; j++)
          Min=fmin(Min,inMin[j]); /* find minimum over a window of length m */
          }
          ptOutMin = *(inMin++);        /* store point comming outMin of the window for future use and move window */
     *(outMin++) = (Min==CSTmin ? NaN : Min);  /* and fill the space with window extreme and move window */

     if(ptOutMax==Max) {        /* if point comining outMax of the window was window's extreme than ... */
     Max=CSTmax;              /* we need to recalculate 'Max' */
     for(j=0; j<m-i-1; j++)
          Max=fmax(Max,inMax[j]); /* find maximum over a window of length m */
     }
     ptOutMax = *(inMax++);        /* store point comming outMax of the window for future use and move window */
     *(outMax++) = (Max==CSTmax ? NaN : Max); /* and fill the space with window extreme and move window */
     }
}

/* the gateway function */
SEXP runminmax(SEXP series, SEXP window)
{
     int w, N;

     // get the length of series
     N = LENGTH(series);

     // window size
     w = asInteger(window);

     if(w > N)
     {
          error("Window cannot be longer than series.");
     }

     // create a C pointer to the output matrix
     SEXP Min = PROTECT(allocVector(REALSXP, N));
     SEXP Max = PROTECT(allocVector(REALSXP, N));

     /*  call the C subroutine */
     runminmax_c(REAL(series), REAL(Min), REAL(Max), N, w);

     // put results in a list
     SEXP list_names = PROTECT(allocVector(STRSXP, 2));
     SET_STRING_ELT(list_names, 0, mkChar("min"));
     SET_STRING_ELT(list_names, 1, mkChar("max"));

     SEXP list = PROTECT(allocVector(VECSXP, 2));
     SET_VECTOR_ELT(list, 0, Min);
     SET_VECTOR_ELT(list, 1, Max);
     setAttrib(list, R_NamesSymbol, list_names);

     UNPROTECT(4);

     return list;
}
