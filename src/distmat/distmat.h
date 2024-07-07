#ifndef DTWCLUST_DISTMAT_HPP_
#define DTWCLUST_DISTMAT_HPP_

#include <memory> // *_ptr

#include <RcppParallel.h>
#include <Rinternals.h>

#include "../utils/utils.h" // id_t

namespace dtwclust {

// -------------------------------------------------------------------------------------------------
/* abstract distmat */
// -------------------------------------------------------------------------------------------------
class Distmat
{
public:
    virtual ~Distmat() {};
    virtual double& operator() (const id_t i, const id_t j) = 0;
    virtual id_t nrow() const = 0;
    virtual id_t ncol() const = 0;
};

// -------------------------------------------------------------------------------------------------
/* concrete factory */
// -------------------------------------------------------------------------------------------------
class DistmatFactory
{
public:
    std::shared_ptr<Distmat> create(const SEXP& MAT_TYPE, const SEXP& D);
};

// -------------------------------------------------------------------------------------------------
/* R matrix distmat (thread-safe) */
// -------------------------------------------------------------------------------------------------
class RDistmat : public Distmat
{
public:
    RDistmat(const SEXP& D);
    double& operator() (const id_t i, const id_t j) override;
    id_t nrow() const override;
    id_t ncol() const override;

private:
    RcppParallel::RMatrix<double> distmat_;
};

} // namespace dtwclust

#endif // DTWCLUST_DISTMAT_HPP_
