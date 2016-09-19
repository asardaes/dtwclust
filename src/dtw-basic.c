// THIS EXPECTS THE SERIES TO SPAN TIME ACROSS ROWS AND DIMENSIONS ACROSS COLUMNS
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <Rdefines.h>
#include "shared.h"

double dtw_no_backtrack(const double *x, const double *y, const int w,
                        const int nx, const int ny, const int dim,
                        const double norm, const double step,
                        double **D)
{
     int i, j;
     int j1, j2;
     double local_cost, global_cost;

     // initialization
     for (i = 0; i <= nx; i++)
     {
          for (j = 0; j <= ny; j++)
               D[i][j] = -1;
     }

     D[1][1] = pow(lnorm(x, y, norm, nx, ny, dim, 0, 0), norm);

     // dynamic programming
     for (i = 1; i <= nx; i++)
     {
          if (w == -1)
          {
               j1 = 1;
               j2 = ny;
          }
          else
          {
               j1 = ceil((double)i * ny / nx - w);
               j2 = floor((double)i * ny / nx + w);

               j1 = j1 > 1 ? j1 : 1;
               j2 = j2 < ny ? j2 : ny;
          }

          for (j = j1; j <= j2; j++)
          {
               if (i == 1 && j == 1) continue;

               local_cost = pow(lnorm(x, y, norm, nx, ny, dim, i-1, j-1), norm);

               if (D[i-1][j] == -1)
                    global_cost = -1;
               else
                    global_cost = D[i-1][j] + local_cost;

               if (D[i][j-1] != -1 && (global_cost == -1 || D[i][j-1] + local_cost < global_cost))
                    global_cost = D[i][j-1] + local_cost;

               if (D[i-1][j-1] != -1 && (global_cost == -1 || D[i-1][j-1] + step*local_cost < global_cost))
                    global_cost = D[i-1][j-1] + step*local_cost;

               D[i][j] = global_cost;
          }
     }

     return pow(D[nx][ny], 1/norm);
}

void backtrack_gcm(double **D, const int nx, const int ny,
                   int *index1, int *index2, int *path)
{
     index1[nx + ny - 1] = nx;
     index2[nx + ny - 1] = ny;
     *path = nx + ny - 2;

     int i = nx;
     int j = ny;
     double triplet[3];
     int temp, which_min;

     double slope = ((double)nx - 1) / (ny - 1);
     double intercept = 1 - slope;

     while (i > 1 || j > 1)
     {
          triplet[0] = (i - 1 > 0) ? D[i-1][j] : DBL_MAX;
          triplet[1] = (i - 1 > 0 && j - 1 > 0) ? D[i-1][j-1] : DBL_MAX;
          triplet[2] = (j - 1 > 0) ? D[i][j-1] : DBL_MAX;

          if (triplet[0] == -1) triplet[0] = DBL_MAX;
          if (triplet[1] == -1) triplet[1] = DBL_MAX;
          if (triplet[2] == -1) triplet[2] = DBL_MAX;

          if (triplet[0] == triplet[1])
          {
               // what leaves me closer to diagonal?
               double diag = slope*j + intercept;

               temp = (i < diag) ? 0 : 1;
          }
          else
          {
               temp = (triplet[0] < triplet[1]) ? 0 : 1;
          }

          if (triplet[temp] == triplet[2])
          {
               // what leaves me closer to diagonal?
               double diag = slope*j + intercept;

               which_min = (i >= diag) ? temp : 2;
          }
          else
          {
               which_min = (triplet[temp] < triplet[2]) ? temp : 2;
          }

          if (which_min < 2 && i > 0) i--;
          if (which_min > 0 && j > 0) j--;

          index1[(*path)] = i;
          index2[(*path)] = j;
          (*path)--;
     }
}

/* the gateway function */
SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP nx, SEXP ny, SEXP dim,
               SEXP norm, SEXP step, SEXP backtrack)
{
     double d;
     int p = 0;
     int m = asInteger(nx);
     int n = asInteger(ny);

     // create gcm
     double **D = (double **)malloc((m + 1) * sizeof(double *));

     for (int i = 0; i <= m; i++)
          D[i] = (double *)malloc((n + 1) * sizeof(double));

     // calculate distance
     d = dtw_no_backtrack(REAL(x), REAL(y), asInteger(window),
                          m, n, asInteger(dim),
                          asReal(norm), asReal(step), D);

     SEXP ret = PROTECT(ScalarReal(d));
     p++;

     if (asLogical(backtrack)) {
          SEXP index1 = PROTECT(allocVector(INTSXP, m + n));
          SEXP index2 = PROTECT(allocVector(INTSXP, m + n));
          p += 2;

          int path = 0;

          backtrack_gcm(D, m, n, INTEGER(index1), INTEGER(index2), &path);

          setAttrib(ret, install("index1"), index1);
          setAttrib(ret, install("index2"), index2);
          setAttrib(ret, install("path"), ScalarInteger(path + 2));
     }

     // free D
     for (int i = 0; i <= m; i++) free(D[i]);
     free(D);

     UNPROTECT(p);
     return ret;
}
