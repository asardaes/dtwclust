#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <Rdefines.h>

int d2s(const int i, const int j, const int nx)
{
     return i + j * (nx + 1);
}

double lnorm(const double *x, const double *y, const double norm,
             const int nx, const int ny, const int dim,
             const int i, const int j)
{
     double res = 0;
     for (int k = 0; k < dim; k++) res += pow(fabs(x[i + nx*k] - y[j + ny*k]), norm);
     res = pow(res, 1/norm);
     return res;
}

double dtw_no_backtrack(const double *x, const double *y, const int w,
                        const int nx, const int ny, const int dim,
                        const double norm, const double step,
                        double *D)
{
     // initialization
     for (int i = 0; i <= nx; i++)
     {
          for (int j = 0; j <= ny; j++)
               D[d2s(i,j,nx)] = -1;
     }

     D[d2s(1,1,nx)] = pow(lnorm(x, y, norm, nx, ny, dim, 0, 0), norm);

     // dynamic programming
     double temp, local_cost, global_cost;

     for (int i = 1; i <= nx; i++)
     {
          int j1, j2;

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

          for (int j = j1; j <= j2; j++)
          {
               if (i == 1 && j == 1) continue;

               local_cost = pow(lnorm(x, y, norm, nx, ny, dim, i-1, j-1), norm);

               temp = D[d2s(i-1,j,nx)];

               if (temp == -1)
                    global_cost = -1;
               else
                    global_cost = temp + local_cost;

               temp = D[d2s(i,j-1,nx)];

               if (temp != -1 && (global_cost == -1 || temp + local_cost < global_cost))
                    global_cost = temp + local_cost;

               temp = D[d2s(i-1,j-1,nx)];

               if (temp != -1 && (global_cost == -1 || temp + step*local_cost < global_cost))
                    global_cost = temp + step*local_cost;

               D[d2s(i,j,nx)] = global_cost;
          }
     }

     return pow(D[(nx+1) * (ny+1) - 1], 1/norm);
}

double dtw_backtrack(const double *x, const double *y, const int w,
                     const int nx, const int ny, const int dim,
                     const double norm, const double step,
                     double *D, int *S, int *index1, int *index2, int *path)
{
     // initialization
     for (int i = 0; i <= nx; i++)
     {
          for (int j = 0; j <= ny; j++)
          {
               D[d2s(i,j,nx)] = -1;
               S[d2s(i,j,nx)] = -1;
          }
     }

     D[d2s(1,1,nx)] = pow(lnorm(x, y, norm, nx, ny, dim, 0, 0), norm);

     // dynamic programming
     double temp, local_cost, global_cost;
     int direction; // 1 = diag, 2 = left, 3 = up

     for (int i = 1; i <= nx; i++)
     {
          int j1, j2;

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

          for (int j = j1; j <= j2; j++)
          {
               if (i == 1 && j == 1) continue;

               local_cost = pow(lnorm(x, y, norm, nx, ny, dim, i-1, j-1), norm);

               temp = D[d2s(i-1,j,nx)];
               direction = 3; // up

               if (temp == -1)
                    global_cost = -1;
               else
                    global_cost = temp + local_cost;

               temp = D[d2s(i,j-1,nx)];

               if (temp != -1 && (global_cost == -1 || temp + local_cost <= global_cost))
               {
                    global_cost = temp + local_cost;
                    direction = 2; // left
               }

               temp = D[d2s(i-1,j-1,nx)];

               if (temp != -1 && (global_cost == -1 || temp + step*local_cost <= global_cost))
               {
                    global_cost = temp + step*local_cost;
                    direction = 1; // diag
               }

               D[d2s(i,j,nx)] = global_cost;
               S[d2s(i,j,nx)] = direction;
          }
     }

     // backtrack
     int i = nx;
     int j = ny;

     index1[0] = nx;
     index2[0] = ny;
     *path = 1;

     while(!(i == 1 && j == 1))
     {
          switch(S[d2s(i,j,nx)])
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
               error("dtw_basic: Invalid direction matrix computed.");
          }

          index1[*path] = i;
          index2[*path] = j;
          (*path)++;
     }

     return pow(D[(nx+1) * (ny+1) - 1], 1/norm);
}

/* the gateway function */
SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP dim,
               SEXP norm, SEXP step, SEXP backtrack,
               SEXP D, SEXP S)
{
     double d;
     int nx = asInteger(m);
     int ny = asInteger(n);

     if (asLogical(backtrack)) {
          SEXP index1 = PROTECT(allocVector(INTSXP, nx + ny));
          SEXP index2 = PROTECT(allocVector(INTSXP, nx + ny));

          int path = 0;

          // calculate distance
          d = dtw_backtrack(REAL(x), REAL(y), asInteger(window),
                            nx, ny, asInteger(dim),
                            asReal(norm), asReal(step),
                            REAL(D), INTEGER(S),
                            INTEGER(index1), INTEGER(index2), &path);

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
          d = dtw_no_backtrack(REAL(x), REAL(y), asInteger(window),
                               nx, ny, asInteger(dim),
                               asReal(norm), asReal(step),
                               REAL(D));

          SEXP ret = PROTECT(ScalarReal(d));

          UNPROTECT(1);
          return ret;
     }
}
