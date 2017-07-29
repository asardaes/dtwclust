#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

// R matrix
void lbi_loop_pairwise(Rcpp::NumericMatrix& dist,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& L, const Rcpp::List& U,
                       Rcpp::NumericVector& L2, Rcpp::NumericVector& U2,
                       Rcpp::NumericVector& H, Rcpp::NumericVector& LB,
                       const int p, const int window, int index)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        Rcpp::NumericVector x(X[i]), y(Y[i]), lower(L[i]), upper(U[i]);
        dist(index++, 0) = lbi_core(x, y, window, p, lower, upper, L2, U2, H, LB);
    }
}

// big.matrix
void lbi_loop_pairwise(MatrixAccessor<double>& dist,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& L, const Rcpp::List& U,
                       Rcpp::NumericVector& L2, Rcpp::NumericVector& U2,
                       Rcpp::NumericVector& H, Rcpp::NumericVector& LB,
                       const int p, const int window, int index)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        Rcpp::NumericVector x(X[i]), y(Y[i]), lower(L[i]), upper(U[i]);
        // bigmemory operator[][] is backwards
        dist[0][index++] = lbi_core(x, y, window, p, lower, upper, L2, U2, H, LB);
    }
}

// =================================================================================================
/* general case */
// =================================================================================================

// R matrix
void lbi_loop_general(Rcpp::NumericMatrix& dist,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& L, const Rcpp::List& U,
                      Rcpp::NumericVector& L2, Rcpp::NumericVector& U2,
                      Rcpp::NumericVector& H, Rcpp::NumericVector& LB,
                      const int p, const int window, int index)
{
    index--;
    for (int j = 0; j < L.length(); j++) {
        Rcpp::NumericVector y(Y[j]), lower(L[j]), upper(U[j]);
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            Rcpp::NumericVector x(X[i]);
            dist(i,index) = lbi_core(x, y, window, p, lower, upper, L2, U2, H, LB);
        }
        index++;
    }
}

// big.matrix
void lbi_loop_general(MatrixAccessor<double>& dist,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& L, const Rcpp::List& U,
                      Rcpp::NumericVector& L2, Rcpp::NumericVector& U2,
                      Rcpp::NumericVector& H, Rcpp::NumericVector& LB,
                      const int p, const int window, int index)
{
    index--;
    for (int j = 0; j < L.length(); j++) {
        Rcpp::NumericVector y(Y[j]), lower(L[j]), upper(U[j]);
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            Rcpp::NumericVector x(X[i]);
            // bigmemory operator[][] is backwards
            dist[index][i] = lbi_core(x, y, window, p, lower, upper, L2, U2, H, LB);
        }
        index++;
    }
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP lbi_loop(SEXP D, SEXP X, SEXP Y, SEXP L, SEXP U,
                         SEXP PAIRWISE, SEXP BIGMAT,
                         SEXP P, SEXP WINDOW, SEXP LEN, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    int len = Rcpp::as<int>(LEN);
    Rcpp::NumericVector L2(len), U2(len), H(len);
    Rcpp::NumericVector LB(len);

    if (Rcpp::as<bool>(BIGMAT)) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> dist(*dist_ptr);

        if (Rcpp::as<bool>(PAIRWISE))
            lbi_loop_pairwise(dist, X, Y, L, U, L2, U2, H, LB,
                              Rcpp::as<int>(P),
                              Rcpp::as<int>(WINDOW),
                              Rcpp::as<int>(ENDPOINTS));
        else
            lbi_loop_general(dist, X, Y, L, U, L2, U2, H, LB,
                             Rcpp::as<int>(P),
                             Rcpp::as<int>(WINDOW),
                             Rcpp::as<int>(ENDPOINTS));

    } else {
        Rcpp::NumericMatrix dist(D);

        if (Rcpp::as<bool>(PAIRWISE))
            lbi_loop_pairwise(dist, X, Y, L, U, L2, U2, H, LB,
                              Rcpp::as<int>(P),
                              Rcpp::as<int>(WINDOW),
                              Rcpp::as<int>(ENDPOINTS));
        else
            lbi_loop_general(dist, X, Y, L, U, L2, U2, H, LB,
                             Rcpp::as<int>(P),
                             Rcpp::as<int>(WINDOW),
                             Rcpp::as<int>(ENDPOINTS));
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
