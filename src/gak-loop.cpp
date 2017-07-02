#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

// R matrix
void gak_loop_pairwise(Rcpp::NumericMatrix& dist,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& distargs,
                       int index, const bool multivariate)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        if (multivariate) {
            Rcpp::NumericMatrix x((SEXP)X[i]);
            Rcpp::NumericMatrix y((SEXP)Y[i]);
            dist(index++, 0) = gak(x, y, distargs);

        } else {
            Rcpp::NumericVector x(X[i]);
            Rcpp::NumericVector y(Y[i]);
            dist(index++, 0) = gak(x, y, distargs);
        }
    }
}

// big.matrix
void gak_loop_pairwise(MatrixAccessor<double>& dist,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& distargs,
                       int index, const bool multivariate)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        if (multivariate) {
            Rcpp::NumericMatrix x((SEXP)X[i]);
            Rcpp::NumericMatrix y((SEXP)Y[i]);
            // bigmemory operator[][] is backwards
            dist[0][index++] = gak(x, y, distargs);

        } else {
            Rcpp::NumericVector x(X[i]);
            Rcpp::NumericVector y(Y[i]);
            // bigmemory operator[][] is backwards
            dist[0][index++] = gak(x, y, distargs);
        }
    }
}

// =================================================================================================
/* symmetric case */
// =================================================================================================

// R matrix
void gak_loop_symmetric(Rcpp::NumericMatrix& dist, const Rcpp::List& X,
                        const Rcpp::List& endpoints, const Rcpp::List& distargs,
                        const bool multivariate)
{
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;

    if (multivariate) {
        while (j < j_end) {
            Rcpp::NumericMatrix y((SEXP)X[j]);

            int i_max;
            if (j == (j_end - 1))
                i_max = i_end;
            else
                i_max = X.length();

            while (i < i_max) {
                R_CheckUserInterrupt();
                Rcpp::NumericMatrix x((SEXP)X[i]);
                double d = gak(x, y, distargs);
                dist(i,j) = d;
                dist(j,i) = d;
                i++;
            }
            j++;
            i = j + 1;
        }

    } else {
        while (j < j_end) {
            Rcpp::NumericVector y(X[j]);

            int i_max;
            if (j == (j_end - 1))
                i_max = i_end;
            else
                i_max = X.length();

            while (i < i_max) {
                R_CheckUserInterrupt();
                Rcpp::NumericVector x(X[i]);
                double d = gak(x, y, distargs);
                dist(i,j) = d;
                dist(j,i) = d;
                i++;
            }
            j++;
            i = j + 1;
        }
    }
}

// big.matrix
void gak_loop_symmetric(MatrixAccessor<double>& dist, const Rcpp::List& X,
                        const Rcpp::List& endpoints, const Rcpp::List& distargs,
                        const bool multivariate)
{
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;

    if (multivariate) {
        while (j < j_end) {
            Rcpp::NumericMatrix y((SEXP)X[j]);

            int i_max;
            if (j == (j_end - 1))
                i_max = i_end;
            else
                i_max = X.length();

            while (i < i_max) {
                R_CheckUserInterrupt();
                Rcpp::NumericMatrix x((SEXP)X[i]);
                double d = gak(x, y, distargs);
                dist[i][j] = d;
                dist[j][i] = d;
                i++;
            }
            j++;
            i = j + 1;
        }

    } else {
        while (j < j_end) {
            Rcpp::NumericVector y(X[j]);

            int i_max;
            if (j == (j_end - 1))
                i_max = i_end;
            else
                i_max = X.length();

            while (i < i_max) {
                R_CheckUserInterrupt();
                Rcpp::NumericVector x(X[i]);
                double d = gak(x, y, distargs);
                dist[i][j] = d;
                dist[j][i] = d;
                i++;
            }
            j++;
            i = j + 1;
        }
    }
}

// =================================================================================================
/* general case */
// =================================================================================================

// R matrix
void gak_loop_general(Rcpp::NumericMatrix& dist, const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& distargs, int index, const bool multivariate)
{
    index--;
    if (multivariate) {
        for (int j = 0; j < Y.length(); j++) {
            Rcpp::NumericMatrix y((SEXP)Y[j]);
            for (int i = 0; i < X.length(); i++) {
                R_CheckUserInterrupt();
                Rcpp::NumericMatrix x((SEXP)X[i]);
                dist(i,index) = gak(x, y, distargs);
            }
            index++;
        }

    } else {
        for (int j = 0; j < Y.length(); j++) {
            Rcpp::NumericVector y(Y[j]);
            for (int i = 0; i < X.length(); i++) {
                R_CheckUserInterrupt();
                Rcpp::NumericVector x(X[i]);
                dist(i,index) = gak(x, y, distargs);
            }
            index++;
        }
    }
}

// big.matrix
void gak_loop_general(MatrixAccessor<double>& dist, const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& distargs, int index, const bool multivariate)
{
    index--;
    if (multivariate) {
        for (int j = 0; j < Y.length(); j++) {
            Rcpp::NumericMatrix y((SEXP)Y[j]);
            for (int i = 0; i < X.length(); i++) {
                R_CheckUserInterrupt();
                Rcpp::NumericMatrix x((SEXP)X[i]);
                // bigmemory operator[][] is backwards
                dist[index][i] = gak(x, y, distargs);
            }
            index++;
        }

    } else {
        for (int j = 0; j < Y.length(); j++) {
            Rcpp::NumericVector y(Y[j]);
            for (int i = 0; i < X.length(); i++) {
                R_CheckUserInterrupt();
                Rcpp::NumericVector x(X[i]);
                // bigmemory operator[][] is backwards
                dist[index][i] = gak(x, y, distargs);
            }
            index++;
        }
    }
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP gak_loop(SEXP D, SEXP X, SEXP Y,
                         SEXP SYMMETRIC, SEXP PAIRWISE,
                         SEXP BIGMAT, SEXP MULTIVARIATE,
                         SEXP DISTARGS, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    if (Rcpp::as<bool>(BIGMAT)) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> dist(*dist_ptr);

        if (Rcpp::as<bool>(PAIRWISE))
            gak_loop_pairwise(dist, X, Y, DISTARGS,
                              Rcpp::as<int>(ENDPOINTS),
                              Rcpp::as<bool>(MULTIVARIATE));
        else if (Rcpp::as<bool>(SYMMETRIC))
            gak_loop_symmetric(dist, X, ENDPOINTS, DISTARGS,
                               Rcpp::as<bool>(MULTIVARIATE));
        else
            gak_loop_general(dist, X, Y, DISTARGS,
                             Rcpp::as<int>(ENDPOINTS),
                             Rcpp::as<bool>(MULTIVARIATE));

    } else {
        Rcpp::NumericMatrix dist(D);

        if (Rcpp::as<bool>(PAIRWISE))
            gak_loop_pairwise(dist, X, Y, DISTARGS,
                              Rcpp::as<int>(ENDPOINTS),
                              Rcpp::as<bool>(MULTIVARIATE));
        else if (Rcpp::as<bool>(SYMMETRIC))
            gak_loop_symmetric(dist, X, ENDPOINTS, DISTARGS,
                               Rcpp::as<bool>(MULTIVARIATE));
        else
            gak_loop_general(dist, X, Y, DISTARGS,
                             Rcpp::as<int>(ENDPOINTS),
                             Rcpp::as<bool>(MULTIVARIATE));
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
