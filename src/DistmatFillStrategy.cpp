#include "dtwclust++.h"
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <memory> // *_ptr

namespace dtwclust {

// =================================================================================================
/* pairwise fill strategy */
// =================================================================================================

void PairwiseDistmatFill::fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                                      const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                      const SEXP& ENDPOINTS, const bool is_bigmat) const
{
    int index = Rcpp::as<int>(ENDPOINTS);
    index--; // R starts at 1, C++ at 0
    if (is_bigmat) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> distmat(*dist_ptr);
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            // bigmemory operator[][] is backwards
            distmat[0][index++] = dist_calculator->calculateDistance(X, Y, i, i);
        }
    }
    else {
        Rcpp::NumericMatrix distmat(D);
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            // lengths of X and Y must be the same, checked in R
            SEXP x = X[i];
            SEXP y = Y[i];
            distmat(index++, 0) = dist_calculator->calculateDistance(X, Y, i, i);
        }
    }
}

// =================================================================================================
/* symmetric fill strategy */
// =================================================================================================

void SymmetricDistmatFill::fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                                       const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                       const SEXP& ENDPOINTS, const bool is_bigmat) const
{
    Rcpp::List endpoints(ENDPOINTS);
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;
    if (is_bigmat) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> distmat(*dist_ptr);
        while (j < j_end) {
            int i_max;
            if (j == (j_end - 1))
                i_max = i_end;
            else
                i_max = X.length();

            while (i < i_max) {
                R_CheckUserInterrupt();
                double d = dist_calculator->calculateDistance(X, X, i, j);
                distmat[i][j] = d;
                distmat[j][i] = d;
                i++;
            }
            j++;
            i = j + 1;
        }
    }
    else {
        Rcpp::NumericMatrix distmat(D);
        while (j < j_end) {
            int i_max;
            if (j == (j_end - 1))
                i_max = i_end;
            else
                i_max = X.length();

            while (i < i_max) {
                R_CheckUserInterrupt();
                double d = dist_calculator->calculateDistance(X, X, i, j);
                distmat(i,j) = d;
                distmat(j,i) = d;
                i++;
            }
            j++;
            i = j + 1;
        }
    }
}

// =================================================================================================
/* general fill strategy */
// =================================================================================================

void GeneralDistmatFill::fillDistmat(const SEXP& D, const Rcpp::List& X, const Rcpp::List& Y,
                                     const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                     const SEXP& ENDPOINTS, const bool is_bigmat) const
{
    int index = Rcpp::as<int>(ENDPOINTS);
    index--; // R starts at 1, C++ at 0
    if (is_bigmat) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> distmat(*dist_ptr);
        for (int j = 0; j < Y.length(); j++) {
            for (int i = 0; i < X.length(); i++) {
                R_CheckUserInterrupt();
                // bigmemory operator[][] is backwards
                distmat[index][i] = dist_calculator->calculateDistance(X, Y, i, j);
            }
            index++;
        }
    }
    else {
        Rcpp::NumericMatrix distmat(D);
        for (int j = 0; j < Y.length(); j++) {
            for (int i = 0; i < X.length(); i++) {
                R_CheckUserInterrupt();
                distmat(i,index) = dist_calculator->calculateDistance(X, Y, i, j);
            }
            index++;
        }
    }
}

} // namespace dtwclust
