#include "fillers.h"

#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppThread.h>

#include "../utils/utils.h" // get_grain, s2d, id_t, interrupt_grain

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
                       const std::shared_ptr<Distmat>& distmat,
                       const int grain)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , interrupt_grain_(interrupt_grain(grain, 10, 1000))
    { }

    // parallel loop across specified range
    void operator()(id_t begin, id_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        for (id_t i = begin; i < end; i++) {
            if (RcppThread::isInterrupted(i % interrupt_grain_ == 0)) break; // nocov
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
    const int interrupt_grain_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void PairwiseFiller::fill() const {
    int size = distmat_->nrow();
    int grain = get_grain(size, num_threads_);
    PairwiseFillWorker fill_worker(dist_calculator_, distmat_, grain);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
    RcppThread::checkUserInterrupt();
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
                      const std::shared_ptr<Distmat>& distmat,
                      const int grain)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , ncols_(distmat->ncol())
        , interrupt_grain_(interrupt_grain(grain, 10, 1000))
    { }

    // parallel loop across specified range
    void operator()(id_t begin, id_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        for (id_t i = begin; i < end; i++) {
            if (RcppThread::isInterrupted()) break; // nocov

            for (id_t j = 0; j < ncols_; j++) {
                if (RcppThread::isInterrupted(j % interrupt_grain_ == 0)) break; // nocov
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
    const id_t ncols_;
    const int interrupt_grain_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void PrimaryFiller::fill() const {
    int size = distmat_->nrow();
    int grain = get_grain(size, num_threads_);
    PrimaryFillWorker fill_worker(dist_calculator_, distmat_, grain);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
    RcppThread::checkUserInterrupt();
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
                        const std::shared_ptr<Distmat>& distmat,
                        const int grain)
        : dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , nrows_(distmat->nrow())
        , interrupt_grain_(interrupt_grain(grain, 10, 1000))
    { }

    // parallel loop across specified range
    void operator()(id_t begin, id_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        id_t i = nrows_;
        id_t j;
        for (id_t id = begin; id < end; id++) {
            if (RcppThread::isInterrupted(id % interrupt_grain_ == 0)) break; // nocov

            if (i >= nrows_ - 1)
                s2d(id, nrows_, i, j); // move to next column
            else
                i++; // same column still valid, only increase row

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
    const id_t nrows_;
    const int interrupt_grain_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void SymmetricFiller::fill() const {
    // number of elements in square matrix without including diagonal
    int nrows = distmat_->nrow();
    int size = nrows * (nrows - 1) / 2;
    int grain = get_grain(size, num_threads_);
    SymmetricFillWorker fill_worker(dist_calculator_, distmat_, grain);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
    RcppThread::checkUserInterrupt();
}

} // namespace dtwclust
