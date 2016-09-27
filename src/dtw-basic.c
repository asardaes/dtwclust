#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <Rdefines.h>

// for cost matrix, in case of window constraint
#define NOT_VISITED -1.0

// for step matrix
#define UP 3.0
#define LEFT 2.0
#define DIAG 1.0

// double to single index, matrices are always vectors in R
int d2s(const int i, const int j, const int nx) __attribute__((always_inline));

int inline d2s(const int i, const int j, const int nx)
{
     return i + j * (nx + 1);
}

// vector norm
double lnorm(const double *x, const double *y, const double norm,
             const int nx, const int ny, const int dim,
             const int i, const int j)
{
     double res = 0;
     for (int k = 0; k < dim; k++) res += pow(fabs(x[i + nx * k] - y[j + ny * k]), norm);
     return (norm == 1) ? res : sqrt(res);
}

// which direction to take in the cost matrix
double which_min(const int i, const int j, const int nx,
                 const double step, const double local_cost,
                 double *D)
{
     volatile double *tuple = (double *)malloc(3 * sizeof(double));

     // DIAG, LEFT, UP
     tuple[0] = (D[d2s(i-1, j-1, nx)] == NOT_VISITED) ? DBL_MAX : D[d2s(i-1, j-1, nx)] + step * local_cost;
     tuple[1] = (D[d2s(i, j-1, nx)] == NOT_VISITED) ? DBL_MAX : D[d2s(i, j-1, nx)] + local_cost;
     tuple[2] = (D[d2s(i-1, j, nx)] == NOT_VISITED) ? DBL_MAX : D[d2s(i-1, j, nx)] + local_cost;

     int min = (tuple[1] < tuple[0]) ? 1 : 0;
     min = (tuple[2] < tuple[min]) ? 2 : min;

     free((double *)tuple);

     return ((double) min + 1.0);
}

// the C code
double dtw_basic_c(const double *x, const double *y, const int w,
                   const int nx, const int ny, const int dim,
                   const double norm, const double step,
                   const int backtrack,
                   double *D, int *index1, int *index2, int *path)
{
     // initialization
     for (int i = 0; i <= nx; i++)
     {
          for (int j = 0; j <= ny; j++)
               D[d2s(i, j, nx)] = NOT_VISITED;
     }

     // first value, must set here to avoid multiplying by step
     D[d2s(1, 1, nx)] = pow(lnorm(x, y, norm, nx, ny, dim, 0, 0), norm);

     // dynamic programming
     double local_cost, global_cost, direction;

     for (int i = 1; i <= nx; i++)
     {
          int j1, j2;

          // adjust limits depending on window
          if (w == -1)
          {
               j1 = 1;
               j2 = ny;
          }
          else
          {
               j1 = ceil((double) i * ny / nx - w);
               j2 = floor((double) i * ny / nx + w);

               j1 = j1 > 1 ? j1 : 1;
               j2 = j2 < ny ? j2 : ny;
          }

          for (int j = j1; j <= j2; j++)
          {
               // very first value already set above
               if (i == 1 && j == 1) continue;

               local_cost = pow(lnorm(x, y, norm, nx, ny, dim, i-1, j-1), norm);

               direction = which_min(i, j, nx, step, local_cost, D);

               if (direction == DIAG)        global_cost = D[d2s(i-1, j-1, nx)] + step * local_cost;
               else if (direction == LEFT)   global_cost = D[d2s(i, j-1, nx)] + local_cost;
               else if (direction == UP)     global_cost = D[d2s(i-1, j, nx)] + local_cost;
               else                          error("dtw_basic: Invalid direction obtained.");

               /*
                * I can use the same matrix to save both cost values and steps taken by shifting
                * the indices left and up for direction. Since the loop advances row-wise, the
                * appropriate values for the cost will be available, and the unnecessary ones are
                * replaced by steps along the way.
                */

               D[d2s(i, j, nx)] = global_cost;
               if (backtrack) D[d2s(i-1, j-1, nx)] = direction;
          }
     }

     // backtrack
     if (backtrack) {
          int i = nx - 1;
          int j = ny - 1;

          // always start at end of series
          index1[0] = nx;
          index2[0] = ny;
          *path = 1;

          while(!(i == 0 && j == 0))
          {
               switch((int) D[d2s(i, j, nx)])
               {
               case 1:
                    i--;
                    j--;
                    break;

               case 2:
                    j--;
                    break;

               case 3:
                    i--;
                    break;

               default:
                    error("dtw_basic: Invalid direction matrix computed. Indices %d and %d.", i+1, j+1);
               }

               index1[*path] = i + 1;
               index2[*path] = j + 1;
               (*path)++;
          }
     }

     return (norm == 1) ? D[(nx+1) * (ny+1) - 1] : sqrt(D[(nx+1) * (ny+1) - 1]);
}

// the gateway function
SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP dim,
               SEXP norm, SEXP step, SEXP backtrack,
               SEXP D)
{
     double d;
     int nx = asInteger(m);
     int ny = asInteger(n);

     if (asLogical(backtrack)) {
          // longest possible path, length will be adjusted in R
          SEXP index1 = PROTECT(allocVector(INTSXP, nx + ny));
          SEXP index2 = PROTECT(allocVector(INTSXP, nx + ny));

          // actual length holder
          int path = 0;

          // calculate distance
          d = dtw_basic_c(REAL(x), REAL(y), asInteger(window),
                          nx, ny, asInteger(dim),
                          asReal(norm), asReal(step), 1,
                          REAL(D), INTEGER(index1), INTEGER(index2), &path);

          // put results in a list
          SEXP list_names = PROTECT(allocVector(STRSXP, 4));
          SET_STRING_ELT(list_names, 0, mkChar("distance"));
          SET_STRING_ELT(list_names, 1, mkChar("index1"));
          SET_STRING_ELT(list_names, 2, mkChar("index2"));
          SET_STRING_ELT(list_names, 3, mkChar("path"));

          SEXP ret = PROTECT(allocVector(VECSXP, 4));
          SET_VECTOR_ELT(ret, 0, PROTECT(ScalarReal(d)));
          SET_VECTOR_ELT(ret, 1, index1);
          SET_VECTOR_ELT(ret, 2, index2);
          SET_VECTOR_ELT(ret, 3, PROTECT(ScalarInteger(path)));
          setAttrib(ret, R_NamesSymbol, list_names);

          UNPROTECT(6);
          return ret;

     } else {
          // calculate distance
          d = dtw_basic_c(REAL(x), REAL(y), asInteger(window),
                          nx, ny, asInteger(dim),
                          asReal(norm), asReal(step), 0,
                          REAL(D), NULL, NULL, NULL);

          SEXP ret = PROTECT(ScalarReal(d));

          UNPROTECT(1);
          return ret;
     }
}
