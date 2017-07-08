#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

// R matrix
void lbk_loop_pairwise(Rcpp::NumericMatrix& dist,
                       const Rcpp::List& X, const Rcpp::List& L, const Rcpp::List& U,
                       Rcpp::NumericVector& H,
                       const int p, int index)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        Rcpp::NumericVector x(X[i]), lower(L[i]), upper(U[i]);
        dist(index++, 0) = lbk_core(x, p, lower, upper, H);
    }
}

// big.matrix
void lbk_loop_pairwise(MatrixAccessor<double>& dist,
                       const Rcpp::List& X, const Rcpp::List& L, const Rcpp::List& U,
                       Rcpp::NumericVector& H,
                       const int p, int index)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        Rcpp::NumericVector x(X[i]), lower(L[i]), upper(U[i]);
        // bigmemory operator[][] is backwards
        dist[0][index++] = lbk_core(x, p, lower, upper, H);
    }
}

// =================================================================================================
/* general case */
// =================================================================================================

// R matrix
void lbk_loop_general(Rcpp::NumericMatrix& dist,
                      const Rcpp::List& X, const Rcpp::List& L, const Rcpp::List& U,
                      Rcpp::NumericVector& H,
                      const int p, int index)
{
    index--;
    for (int j = 0; j < L.length(); j++) {
        Rcpp::NumericVector lower(L[j]), upper(U[j]);
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            Rcpp::NumericVector x(X[i]);
            dist(i,index) = lbk_core(x, p, lower, upper, H);
        }
        index++;
    }
}

// big.matrix
void lbk_loop_general(MatrixAccessor<double>& dist,
                      const Rcpp::List& X, const Rcpp::List& L, const Rcpp::List& U,
                      Rcpp::NumericVector& H,
                      const int p, int index)
{
    index--;
    for (int j = 0; j < L.length(); j++) {
        Rcpp::NumericVector lower(L[j]), upper(U[j]);
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            Rcpp::NumericVector x(X[i]);
            // bigmemory operator[][] is backwards
            dist[index][i] = lbk_core(x, p, lower, upper, H);
        }
        index++;
    }
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP lbk_loop(SEXP D, SEXP X, SEXP L, SEXP U,
                         SEXP PAIRWISE, SEXP BIGMAT,
                         SEXP P, SEXP LEN, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    Rcpp::NumericVector H(Rcpp::as<int>(LEN));

    if (Rcpp::as<bool>(BIGMAT)) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> dist(*dist_ptr);

        if (Rcpp::as<bool>(PAIRWISE))
            lbk_loop_pairwise(dist, X, L, U, H,
                              Rcpp::as<int>(P),
                              Rcpp::as<int>(ENDPOINTS));
        else
            lbk_loop_general(dist, X, L, U, H,
                             Rcpp::as<int>(P),
                             Rcpp::as<int>(ENDPOINTS));

    } else {
        Rcpp::NumericMatrix dist(D);

        if (Rcpp::as<bool>(PAIRWISE))
            lbk_loop_pairwise(dist, X, L, U, H,
                              Rcpp::as<int>(P),
                              Rcpp::as<int>(ENDPOINTS));
        else
            lbk_loop_general(dist, X, L, U, H,
                             Rcpp::as<int>(P),
                             Rcpp::as<int>(ENDPOINTS));
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
