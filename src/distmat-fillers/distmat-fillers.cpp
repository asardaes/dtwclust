#include "distmat-fillers.h"

#include <cstddef> // std::size_t
#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../utils/utils.h" // get_grain, s2d

namespace dtwclust {

// =================================================================================================
/* Factory methods */
// =================================================================================================

std::shared_ptr<DistmatFiller>
DistmatFillerFactory::create(const SEXP& FILL_TYPE,
                             const SEXP& NUM_THREADS,
                             std::shared_ptr<Distmat>& distmat,
                             const std::shared_ptr<DistanceCalculator>& dist_calculator)
{
    std::string type = Rcpp::as<std::string>(FILL_TYPE);
    if (type == "PAIRWISE") {
        return std::make_shared<PairwiseFiller>(distmat, dist_calculator, NUM_THREADS);
    }
    else if (type == "SYMMETRIC") {
        return std::make_shared<SymmetricFiller>(distmat, dist_calculator, NUM_THREADS);
    }
    else {
        return std::make_shared<PrimaryFiller>(distmat, dist_calculator, NUM_THREADS);
    }
}

// =================================================================================================
/* multi-threaded pairwise filler
 *   This filler assumes a column vector has been received in *distmat_, and that the distance
 *   calculator has the same amount of series in X and Y.
 */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
PairwiseFiller::PairwiseFiller(std::shared_ptr<Distmat>& distmat,
                               const std::shared_ptr<DistanceCalculator>& dist_calculator,
                               const SEXP& NUM_THREADS)
    : DistmatFiller(distmat, dist_calculator, NUM_THREADS)
{ }

// -------------------------------------------------------------------------------------------------
/* parallel worker */
// -------------------------------------------------------------------------------------------------
class PairwiseFillWorker : public RcppParallel::Worker {
public:
    // constructor
    PairwiseFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                       const std::shared_ptr<Distmat>& distmat)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();
        // fill distances
        for (std::size_t i = begin; i < end; i++) {
            (*distmat_)(i,0) = dist_calculator->calculate(i,i);
        }
        // delete local calculator
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    std::shared_ptr<Distmat> distmat_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void PairwiseFiller::fill() const {
    PairwiseFillWorker fill_worker(dist_calculator_, distmat_);
    int size = distmat_->nrow();
    int grain = get_grain(size, num_threads_);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
}

// =================================================================================================
/* multi-threaded primary filler
 *   This filler makes no assumptions about the dimensions of *distmat_, but makes the parallel
 *   chunks based on the number of rows.
 */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
PrimaryFiller::PrimaryFiller(std::shared_ptr<Distmat>& distmat,
                             const std::shared_ptr<DistanceCalculator>& dist_calculator,
                             const SEXP& NUM_THREADS)
    : DistmatFiller(distmat, dist_calculator, NUM_THREADS)
{ }

// -------------------------------------------------------------------------------------------------
/* parallel worker */
// -------------------------------------------------------------------------------------------------
class PrimaryFillWorker : public RcppParallel::Worker {
public:
    // constructor
    PrimaryFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                      const std::shared_ptr<Distmat>& distmat)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , ncols_(distmat->ncol())
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();
        // fill distances
        for (std::size_t i = begin; i < end; i++) {
            for (int j = 0; j < ncols_; j++) {
                (*distmat_)(i,j) = dist_calculator->calculate(i,j);
            }
        }
        // delete local calculator
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    std::shared_ptr<Distmat> distmat_;
    int ncols_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void PrimaryFiller::fill() const {
    PrimaryFillWorker fill_worker(dist_calculator_, distmat_);
    int size = distmat_->nrow();
    int grain = get_grain(size, num_threads_);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
}

// =================================================================================================
/* multi-threaded symmetric filler
 *   This filler assumes a square matrix has been received in *distmat_
 */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
SymmetricFiller::SymmetricFiller(std::shared_ptr<Distmat>& distmat,
                                 const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                 const SEXP& NUM_THREADS)
    : DistmatFiller(distmat, dist_calculator, NUM_THREADS)
{ }

// -------------------------------------------------------------------------------------------------
/* parallel worker */
// -------------------------------------------------------------------------------------------------
class SymmetricFillWorker : public RcppParallel::Worker {
public:
    // constructor
    SymmetricFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                        const std::shared_ptr<Distmat>& distmat)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , nrows_(distmat->nrow())
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();
        // fill distances
        int i, j;
        for (std::size_t id = begin; id < end; id++) {
            s2d(id, nrows_, i, j);
            double dist = dist_calculator->calculate(i,j);
            (*distmat_)(i,j) = dist;
            (*distmat_)(j,i) = dist;
        }
        // delete local calculator
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    std::shared_ptr<Distmat> distmat_;
    int nrows_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void SymmetricFiller::fill() const {
    SymmetricFillWorker fill_worker(dist_calculator_, distmat_);
    // number of elements in square matrix without including diagonal
    int nrows = distmat_->nrow();
    int size = nrows * (nrows + 1) / 2 - nrows;
    int grain = get_grain(size, num_threads_);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
}

} // namespace dtwclust
