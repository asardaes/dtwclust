#ifndef DTWCLUST_TSTSLIST_HPP_
#define DTWCLUST_TSTSLIST_HPP_

#include <complex>
#include <memory> // make_shared, shared_ptr
#include <utility> // move
#include <vector>

#include <RcppArmadillo.h>

namespace dtwclust {

// =================================================================================================
/** Thread-Safe Time-Series List
 *
 *  It must be constructed outside the thread though.
 *
 *  TSTSType must be either arma::mat or arma::cx_mat.
 */

// primary
template<typename TSTSType>
class TSTSList
{
public:
    // constructors
    TSTSList() {} // dummy
    TSTSList(const Rcpp::List& series)
        : series_(std::make_shared<std::vector<TSTSType>>())
    {
        for (const SEXP& x : series) {
            if (Rf_isMatrix(x)) {
                Rcpp::NumericMatrix x_rcpp(x);
                TSTSType x_arma(&x_rcpp[0], x_rcpp.nrow(), x_rcpp.ncol(), false, true);
                series_->push_back(std::move(x_arma));
            }
            else {
                Rcpp::NumericVector x_rcpp(x);
                TSTSType x_arma(&x_rcpp[0], x_rcpp.length(), 1, false, true);
                series_->push_back(std::move(x_arma));
            }
        }
    }
    // operator[]
    TSTSType& operator[](const int i) { return (*series_)[i]; }
    const TSTSType& operator[](const int i) const { return (*series_)[i]; }
    // length
    int length() const { return series_->size(); }
    // iterators
    typename std::vector<TSTSType>::iterator begin() { return series_->begin(); }
    typename std::vector<TSTSType>::iterator end() { return series_->end(); }
    typename std::vector<TSTSType>::const_iterator begin() const { return series_->cbegin(); }
    typename std::vector<TSTSType>::const_iterator end() const { return series_->cend(); }

private:
    std::shared_ptr<std::vector<TSTSType>> series_;
};

// specialization for complex case
template<>
class TSTSList<arma::cx_mat> {
public:
    // constructors
    TSTSList() {} // dummy
    TSTSList(const Rcpp::List& series)
        : series_(std::make_shared<std::vector<arma::cx_mat>>())
    {
        // see http://rcpp-devel.r-forge.r-project.narkive.com/o5ubHVos/multiplication-of-complexvector
        for (const SEXP& x : series) {
            // only complex vectors expected right now
            Rcpp::ComplexVector x_rcpp(x);
            arma::cx_mat x_arma(reinterpret_cast<std::complex<double>*>(&x_rcpp[0]),
                                x_rcpp.length(), 1, false, true);
            series_->push_back(std::move(x_arma));
        }
    }
    // operator[]
    arma::cx_mat& operator[](const int i) { return (*series_)[i]; }
    const arma::cx_mat& operator[](const int i) const { return (*series_)[i]; }
    // length
    int length() const { return series_->size(); }
    // iterators
    typename std::vector<arma::cx_mat>::iterator begin() { return series_->begin(); }
    typename std::vector<arma::cx_mat>::iterator end() { return series_->end(); }
    typename std::vector<arma::cx_mat>::const_iterator begin() const { return series_->cbegin(); }
    typename std::vector<arma::cx_mat>::const_iterator end() const { return series_->cend(); }

private:
    std::shared_ptr<std::vector<arma::cx_mat>> series_;
};

} // namespace dtwclust

#endif // DTWCLUST_TSTSLIST_HPP_
