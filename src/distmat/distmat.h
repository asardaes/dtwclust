#ifndef DTWCLUST_DISTMAT_HPP_
#define DTWCLUST_DISTMAT_HPP_

#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <memory> // *_ptr

namespace dtwclust {

// =================================================================================================
/* Distmat (base + factory + concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distmat */
// -------------------------------------------------------------------------------------------------
class Distmat
{
public:
    virtual ~Distmat() {};
    virtual double& operator() (const int i, const int j) = 0;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistmatFactory
{
public:
    std::shared_ptr<Distmat>
    create(const SEXP& MAT_TYPE, const SEXP& D);
};

// -------------------------------------------------------------------------------------------------
/* R matrix distmat */
// -------------------------------------------------------------------------------------------------
class RDistmat : public Distmat
{
public:
    RDistmat(const SEXP& D);
    double& operator() (const int i, const int j) override;

private:
    Rcpp::NumericMatrix distmat_;
};

// -------------------------------------------------------------------------------------------------
/* bigmemory big.matrix distmat */
// -------------------------------------------------------------------------------------------------
class BigmemoryDistmat : public Distmat
{
public:
    BigmemoryDistmat(const SEXP& D);
    double& operator() (const int i, const int j) override;

private:
    MatrixAccessor<double> distmat_;
};

} // namespace dtwclust

#endif // DTWCLUST_DISTMAT_HPP_
