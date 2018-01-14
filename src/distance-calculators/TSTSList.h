#ifndef DTWCLUST_TSTSLIST_HPP_
#define DTWCLUST_TSTSLIST_HPP_

#include <complex>
#include <type_traits> // conditional, is_same
#include <vector>

#include <RcppArmadillo.h> // arma:: referenced here
#include <RcppParallel.h>

namespace dtwclust {

// =================================================================================================
/* Thread-Safe Time-Series List
 *   It should be constructed outside the thread though.
 */
// =================================================================================================

// primary one
template<typename SeriesType>
class TSTSList
{
public:
    typedef typename std::conditional<
        std::is_same<SeriesType, Rcpp::NumericVector>::value,
        RcppParallel::RVector<double>,
        RcppParallel::RMatrix<double>
    >::type TSTSType;
    // constructors
    TSTSList() {} // dummy
    TSTSList(const Rcpp::List& series)
    {
        for (const SEXP& x : series) {
            SeriesType x_rcpp(x);
            series_.push_back(TSTSType(x_rcpp));
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return series_[i]; }
    const TSTSType& operator[](const int i) const { return series_[i]; }

private:
    std::vector<TSTSType> series_;
};

// specialization for complex vectors
template<>
class TSTSList<Rcpp::ComplexVector>
{
public:
    typedef arma::cx_vec TSTSType;
    // constructors
    TSTSList() {} // dummy
    TSTSList(const Rcpp::List& series)
    {
        for (const SEXP& x : series) {
            series_.push_back(arma::cx_vec(Rcpp::ComplexVector(x)));
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return series_[i]; }
    const TSTSType& operator[](const int i) const { return series_[i]; }

private:
    std::vector<TSTSType> series_;
};

} // namespace dtwclust

#endif // DTWCLUST_TSTSLIST_HPP_
