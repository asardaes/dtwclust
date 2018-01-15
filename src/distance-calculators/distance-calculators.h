#ifndef DTWCLUST_DISTANCE_CALCULATORS_HPP_
#define DTWCLUST_DISTANCE_CALCULATORS_HPP_

#include <memory> // *_ptr
#include <string>

#include <RcppArmadillo.h>

namespace dtwclust {

// =================================================================================================
/* DistanceCalculator (base + factory) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distance calculator */
// -------------------------------------------------------------------------------------------------
class DistanceCalculator
{
public:
    virtual ~DistanceCalculator() {}
    virtual double calculate(const int i, const int j) = 0;
    // a clone method to make life easier when copying objects in each thread
    virtual DistanceCalculator* clone() const = 0;
    // helpers for distmat filler
    virtual int xLimit() const = 0;
    virtual int yLimit() const = 0;

protected:
    int maxLength(const Rcpp::List& list, const bool is_multivariate) const {
        int max_len = 0;
        for (const SEXP& series : list) {
            if (is_multivariate) {
                Rcpp::NumericMatrix x(series);
                int this_len = x.nrow();
                if (this_len > max_len) max_len = this_len;
            }
            else {
                Rcpp::NumericVector x(series);
                int this_len = x.length();
                if (this_len > max_len) max_len = this_len;
            }
        }
        return max_len;
    }
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistanceCalculatorFactory
{
public:
    std::shared_ptr<DistanceCalculator> create(const SEXP& DIST, const SEXP& DIST_ARGS,
                                               const SEXP& X, const SEXP& Y);
    std::shared_ptr<DistanceCalculator> create(const std::string& dist, const SEXP& DIST_ARGS,
                                               const SEXP& X, const SEXP& Y);
};

} // namespace dtwclust

#endif // DTWCLUST_DISTANCE_CALCULATORS_HPP_
