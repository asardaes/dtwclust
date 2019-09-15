#ifndef DTWCLUST_DISTANCE_CALCULATORS_HPP_
#define DTWCLUST_DISTANCE_CALCULATORS_HPP_

#include <memory> // *_ptr
#include <string>

#include <RcppArmadillo.h> // arma:: referenced here

#include "../utils/SurrogateMatrix.h"
#include "../utils/TSTSList.h"
#include "../utils/utils.h" // id_t

namespace dtwclust {

// =================================================================================================
/* DistanceCalculator (abstract + factory + concretes) */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* abstract distance calculator */
// -------------------------------------------------------------------------------------------------
class DistanceCalculator
{
public:
    virtual ~DistanceCalculator() {}
    virtual double calculate(const id_t i, const id_t j) = 0;
    // a clone method to make life easier when copying objects in each thread
    virtual DistanceCalculator* clone() const = 0;

protected:
    int maxLength(const TSTSList<arma::mat>& list) const {
        unsigned int max_len = 0;
        for (const arma::mat& x : list) {
            if (x.n_rows > max_len) max_len = x.n_rows;
        }
        return static_cast<int>(max_len);
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

// -------------------------------------------------------------------------------------------------
/* dtw_basic calculator */
// -------------------------------------------------------------------------------------------------
class DtwBasicCalculator : public DistanceCalculator
{
public:
    DtwBasicCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const id_t i, const id_t j) override;
    DtwBasicCalculator* clone() const override;

private:
    // method calculate
    double calculate(const arma::mat& x, const arma::mat& y);
    // input parameters
    int window_;
    double norm_, step_;
    bool normalize_;
    bool sqrt_dist_;
    // input series
    TSTSList<arma::mat> x_, y_;
    // helper "matrix"
    SurrogateMatrix<double> lcm_;
    // to dimension gcm_
    int max_len_y_;
};

// -------------------------------------------------------------------------------------------------
/* gak calculator */
// -------------------------------------------------------------------------------------------------
class GakCalculator : public DistanceCalculator
{
public:
    GakCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const id_t i, const id_t j) override;
    GakCalculator* clone() const override;

private:
    // method calculate
    double calculate(const arma::mat& x, const arma::mat& y);
    // input parameters
    double sigma_;
    int window_;
    // input series
    TSTSList<arma::mat> x_, y_;
    // helper "matrix"
    SurrogateMatrix<double> logs_;
    // to dimension logs_
    int max_len_x_, max_len_y_;
};

// -------------------------------------------------------------------------------------------------
/* lb_improved calculator */
// -------------------------------------------------------------------------------------------------
class LbiCalculator : public DistanceCalculator
{
public:
    LbiCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const id_t i, const id_t j) override;
    LbiCalculator* clone() const override;

private:
    double calculate(const arma::mat& x, const arma::mat& y,
                     const arma::mat& lower_envelope, const arma::mat& upper_envelope);
    int p_, len_;
    unsigned int window_;
    TSTSList<arma::mat> x_, y_, lower_envelopes_, upper_envelopes_;
    SurrogateMatrix<double> H_, L2_, U2_, LB_;
};

// -------------------------------------------------------------------------------------------------
/* lb_keogh calculator */
// -------------------------------------------------------------------------------------------------
class LbkCalculator : public DistanceCalculator
{
public:
    LbkCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const id_t i, const id_t j) override;
    LbkCalculator* clone() const override;

private:
    double calculate(const arma::mat& x,
                     const arma::mat& lower_envelope, const arma::mat& upper_envelope);
    int p_, len_;
    TSTSList<arma::mat> x_, lower_envelopes_, upper_envelopes_;
    SurrogateMatrix<double> H_;
};

// -------------------------------------------------------------------------------------------------
/* sbd calculator */
// -------------------------------------------------------------------------------------------------
class SbdCalculator : public DistanceCalculator
{
public:
    SbdCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const id_t i, const id_t j) override;
    SbdCalculator* clone() const override;

private:
    double calculate(const arma::mat& x, const arma::mat& y,
                     const arma::cx_mat& fftx, const arma::cx_mat& ffty);
    TSTSList<arma::mat> x_, y_;
    TSTSList<arma::cx_mat> fftx_, ffty_;
    arma::vec cc_seq_truncated_;
    int fftlen_;
};

// -------------------------------------------------------------------------------------------------
/* soft-DTW calculator */
// -------------------------------------------------------------------------------------------------
class SdtwCalculator : public DistanceCalculator
{
public:
    SdtwCalculator(const SEXP& DIST_ARGS, const SEXP& X, const SEXP& Y);
    double calculate(const id_t i, const id_t j) override;
    SdtwCalculator* clone() const override;

private:
    // method calculate
    double calculate(const arma::mat& x, const arma::mat& y);
    // input parameters
    double gamma_;
    // input series
    TSTSList<arma::mat> x_, y_;
    // helper "matrix"
    SurrogateMatrix<double> cm_;
    // to dimension cm_
    int max_len_x_, max_len_y_;
};

} // namespace dtwclust

#endif // DTWCLUST_DISTANCE_CALCULATORS_HPP_
