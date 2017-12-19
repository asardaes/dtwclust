#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

// R matrix
void sdtw_loop_pairwise(Rcpp::NumericMatrix& dist,
                        const Rcpp::List& X, const Rcpp::List& Y,
                        const SEXP& GAMMA, SEXP& COSTMAT,
                        int index, const SEXP& MV)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        SEXP x = X[i];
        SEXP y = Y[i];
        double d = Rcpp::as<double>(soft_dtw(x, y, GAMMA, COSTMAT, R_NilValue, MV));
        dist(index++, 0) = d;
    }
}

// big.matrix
void sdtw_loop_pairwise(MatrixAccessor<double>& dist,
                        const Rcpp::List& X, const Rcpp::List& Y,
                        const SEXP& GAMMA, SEXP& COSTMAT,
                        int index, const SEXP& MV)
{
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        SEXP x = X[i];
        SEXP y = Y[i];
        double d = Rcpp::as<double>(soft_dtw(x, y, GAMMA, COSTMAT, R_NilValue, MV));
        // bigmemory operator[][] is backwards
        dist[0][index++] = d;
    }
}

// =================================================================================================
/* symmetric case */
// =================================================================================================

// R matrix
void sdtw_loop_symmetric(Rcpp::NumericMatrix& dist, const Rcpp::List& X,
                         const SEXP& GAMMA, SEXP& COSTMAT,
                         const Rcpp::List& endpoints, const SEXP& MV)
{
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];
    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        SEXP y = X[j];

        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = X.length();

        while (i < i_max) {
            R_CheckUserInterrupt();
            SEXP x = X[i];
            double d = Rcpp::as<double>(soft_dtw(x, y, GAMMA, COSTMAT, R_NilValue, MV));
            dist(i,j) = d;
            dist(j,i) = d;
            i++;
        }
        j++;
        i = j + 1;
    }
}

// big.matrix
void sdtw_loop_symmetric(MatrixAccessor<double>& dist, const Rcpp::List& X,
                         const SEXP& GAMMA, SEXP& COSTMAT,
                         const Rcpp::List& endpoints, const SEXP& MV)
{
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];
    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        SEXP y = X[j];

        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = X.length();

        while (i < i_max) {
            R_CheckUserInterrupt();
            SEXP x = X[i];
            double d = Rcpp::as<double>(soft_dtw(x, y, GAMMA, COSTMAT, R_NilValue, MV));
            dist[i][j] = d;
            dist[j][i] = d;
            i++;
        }
        j++;
        i = j + 1;
    }
}

// =================================================================================================
/* general case */
// =================================================================================================

// R matrix
void sdtw_loop_general(Rcpp::NumericMatrix& dist, const Rcpp::List& X, const Rcpp::List& Y,
                       const SEXP& GAMMA, SEXP& COSTMAT, int index, const SEXP& MV)
{
    index--;
    for (int j = 0; j < Y.length(); j++) {
        SEXP y = Y[j];
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            SEXP x = X[i];
            double d = Rcpp::as<double>(soft_dtw(x, y, GAMMA, COSTMAT, R_NilValue, MV));
            dist(i,index) = d;
        }
        index++;
    }
}

// big.matrix
void sdtw_loop_general(MatrixAccessor<double>& dist, const Rcpp::List& X, const Rcpp::List& Y,
                       const SEXP& GAMMA, SEXP& COSTMAT, int index, const SEXP& MV)
{
    index--;
    for (int j = 0; j < Y.length(); j++) {
        SEXP y = Y[j];
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            SEXP x = X[i];
            double d = Rcpp::as<double>(soft_dtw(x, y, GAMMA, COSTMAT, R_NilValue, MV));
            // bigmemory operator[][] is backwards
            dist[index][i] = d;
        }
        index++;
    }
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP sdtw_loop(SEXP D, SEXP X, SEXP Y, SEXP DISTARGS,
                          SEXP SYMMETRIC, SEXP PAIRWISE,
                          SEXP BIGMAT, SEXP MV,
                          SEXP ENDPOINTS)
{
    BEGIN_RCPP
    Rcpp::List dots(DISTARGS);
    SEXP GAMMA = dots["gamma"];
    SEXP COSTMAT = dots["cm"];

    if (Rcpp::as<bool>(BIGMAT)) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> dist(*dist_ptr);

        if (Rcpp::as<bool>(PAIRWISE))
            sdtw_loop_pairwise(dist, X, Y, GAMMA, COSTMAT,
                               Rcpp::as<int>(ENDPOINTS),
                               MV);
        else if (Rcpp::as<bool>(SYMMETRIC))
            sdtw_loop_symmetric(dist, X, GAMMA, COSTMAT,
                                ENDPOINTS,
                                MV);
        else
            sdtw_loop_general(dist, X, Y, GAMMA, COSTMAT,
                              Rcpp::as<int>(ENDPOINTS),
                              MV);

    } else {
        Rcpp::NumericMatrix dist(D);

        if (Rcpp::as<bool>(PAIRWISE))
            sdtw_loop_pairwise(dist, X, Y, GAMMA, COSTMAT,
                               Rcpp::as<int>(ENDPOINTS),
                               MV);
        else if (Rcpp::as<bool>(SYMMETRIC))
            sdtw_loop_symmetric(dist, X, GAMMA, COSTMAT,
                                ENDPOINTS,
                                MV);
        else
            sdtw_loop_general(dist, X, Y, GAMMA, COSTMAT,
                              Rcpp::as<int>(ENDPOINTS),
                              MV);
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
