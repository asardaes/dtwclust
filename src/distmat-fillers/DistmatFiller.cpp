#include "distmat-fillers.h"

#include <RcppArmadillo.h>

#include "../distance-calculators/distance-calculators.h"
#include "../distmat/distmat.h"

namespace dtwclust {

// =================================================================================================
/* pairwise */
// =================================================================================================

void PairwiseDistmatFiller::fill() const
{
    int index = Rcpp::as<int>(endpoints_);
    index--; // R starts at 1, C++ at 0
    for (int i = 0; i < dist_calculator_->xLimit(); i++) {
        if (i % 1000 == 0) Rcpp::checkUserInterrupt();
        (*distmat_)(index++, 0) = dist_calculator_->calculate(i,i);
    }
}

// =================================================================================================
/* symmetric */
// =================================================================================================

void SymmetricDistmatFiller::fill() const
{
    Rcpp::List endpoints(endpoints_);
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        Rcpp::checkUserInterrupt();

        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = dist_calculator_->xLimit();

        while (i < i_max) {
            double d = dist_calculator_->calculate(i,j);
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

void GeneralDistmatFiller::fill() const
{
    int index = Rcpp::as<int>(endpoints_);
    index--; // R starts at 1, C++ at 0
    for (int j = 0; j < dist_calculator_->yLimit(); j++) {
        Rcpp::checkUserInterrupt();
        for (int i = 0; i < dist_calculator_->xLimit(); i++) {
            (*distmat_)(i, index) = dist_calculator_->calculate(i,j);
        }
        index++;
    }
}

} // namespace dtwclust
