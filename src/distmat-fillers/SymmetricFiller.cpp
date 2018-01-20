#include "concrete-fillers.h"

#include <cstddef>
#include <memory>

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../utils/utils.h" // get_grain, s2d

namespace dtwclust {

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
        , nrows_(dist_calculator_->xLimit())
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
void SymmetricFiller::fill() const
{
    SymmetricFillWorker fill_worker(dist_calculator_, distmat_);
    // number of elements in square matrix without including diagonal
    int nrows = dist_calculator_->xLimit();
    int size = nrows * (nrows + 1) / 2 - nrows;
    int grain = get_grain(size, num_threads_);
    RcppParallel::parallelFor(0, size, fill_worker, grain);
}

} // namespace dtwclust
