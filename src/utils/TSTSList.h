#ifndef DTWCLUST_TSTSLIST_HPP_
#define DTWCLUST_TSTSLIST_HPP_

#include <complex>
#include <memory> // make_shared, shared_ptr
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
        : series_(std::make_shared<std::vector<TSTSType>>())
    {
        for (const SEXP& x : series) {
            SeriesType x_rcpp(x);
            series_->push_back(TSTSType(x_rcpp));
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return (*series_)[i]; }
    const TSTSType& operator[](const int i) const { return (*series_)[i]; }
    // length
    int length() const { return series_->size(); }

private:
    std::shared_ptr<std::vector<TSTSType>> series_;
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
        : series_(std::make_shared<std::vector<TSTSType>>())
    {
        for (const SEXP& x : series) {
            Rcpp::ComplexVector x_rcpp(x);
            // see http://rcpp-devel.r-forge.r-project.narkive.com/o5ubHVos/multiplication-of-complexvector
            arma::cx_vec x_arma(reinterpret_cast<std::complex<double>*>(x_rcpp.begin()),
                                x_rcpp.length(), false, true);
            series_->push_back(x_arma);
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return (*series_)[i]; }
    const TSTSType& operator[](const int i) const { return (*series_)[i]; }
    // length
    int length() const { return series_->size(); }

private:
    std::shared_ptr<std::vector<TSTSType>> series_;
};

} // namespace dtwclust

#endif // DTWCLUST_TSTSLIST_HPP_
