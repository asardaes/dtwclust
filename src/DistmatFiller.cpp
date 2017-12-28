#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* pairwise */
// =================================================================================================

void PairwiseDistmatFiller::fill(const Rcpp::List& X, const Rcpp::List& Y) const
{
    int index = Rcpp::as<int>(endpoints_);
    index--; // R starts at 1, C++ at 0
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        (*distmat_)(index++, 0) = dist_calculator_->calculate(X, Y, i, i);
    }
}

// =================================================================================================
/* symmetric */
// =================================================================================================

void SymmetricDistmatFiller::fill(const Rcpp::List& X, const Rcpp::List& Y) const
{
    Rcpp::List endpoints(endpoints_);
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = X.length();

        while (i < i_max) {
            R_CheckUserInterrupt();
            double d = dist_calculator_->calculate(X, X, i, j);
            (*distmat_)(i,j) = d;
            (*distmat_)(j,i) = d;
            i++;
        }
        j++;
        i = j + 1;
    }
}

// =================================================================================================
/* general */
// =================================================================================================

void GeneralDistmatFiller::fill(const Rcpp::List& X, const Rcpp::List& Y) const
{
    int index = Rcpp::as<int>(endpoints_);
    index--; // R starts at 1, C++ at 0
    for (int j = 0; j < Y.length(); j++) {
        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();
            (*distmat_)(i, index) = dist_calculator_->calculate(X, Y, i, j);
        }
        index++;
    }
}

} // namespace dtwclust
