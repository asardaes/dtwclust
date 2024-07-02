#include "fillers.h"

#include <cmath> // std::sqrt
#include <memory> // make_shared
#include <string>

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../utils/ParallelWorker.h"
#include "../utils/utils.h" // get_grain, s2d, id_t

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
    // else if (type == "SYMMETRIC") {
    //     return std::make_shared<SymmetricFiller>(distmat, dist_calculator, NUM_THREADS);
    // }
    else if (type == "LOWER_TRIANGULAR") {
        return std::make_shared<LowerTriangularFiller>(distmat, dist_calculator, NUM_THREADS);
    }
    else if (type == "LOWER_TRIANGULAR_DIAGONAL") {
        return std::make_shared<LowerTriangularDiagonalFiller>(distmat, dist_calculator, NUM_THREADS);
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
class PairwiseFillWorker : public ParallelWorker {
public:
    // constructor
    PairwiseFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                       const std::shared_ptr<Distmat>& distmat,
                       const int grain)
        : ParallelWorker(grain, 10, 1000)
        , dist_calculator_(dist_calculator)
        , distmat_(distmat)
    { }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        for (id_t i = begin; i < end; i++) {
            if (is_interrupted(i)) break; // nocov
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
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void PairwiseFiller::fill() const {
    int size = distmat_->nrow();
    int grain = get_grain(size, num_threads_);
    PairwiseFillWorker fill_worker(dist_calculator_, distmat_, grain);
    parallel_for(0, size, fill_worker, grain);
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
class PrimaryFillWorker : public ParallelWorker {
public:
    // constructor
    PrimaryFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                      const std::shared_ptr<Distmat>& distmat,
                      const int grain)
        : ParallelWorker(grain, 10, 1000)
        , dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , ncols_(distmat->ncol())
    { }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        for (id_t i = begin; i < end; i++) {
            if (is_interrupted()) break; // nocov

            for (id_t j = 0; j < ncols_; j++) {
                if (is_interrupted(j)) break; // nocov
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
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void PrimaryFiller::fill() const {
    int size = distmat_->nrow();
    int grain = get_grain(size, num_threads_);
    PrimaryFillWorker fill_worker(dist_calculator_, distmat_, grain);
    parallel_for(0, size, fill_worker, grain);
}

// =================================================================================================
/* multi-threaded symmetric filler
 *   This filler assumes a square matrix has been received in *distmat_
 */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
// SymmetricFiller::SymmetricFiller(std::shared_ptr<Distmat>& distmat,
//                                  const std::shared_ptr<DistanceCalculator>& dist_calculator,
//                                  const SEXP& NUM_THREADS)
//     : DistmatFiller(distmat, dist_calculator, NUM_THREADS)
// { }

// -------------------------------------------------------------------------------------------------
/* parallel worker */
// -------------------------------------------------------------------------------------------------
class SymmetricFillWorker : public ParallelWorker {
public:
    // constructor
    SymmetricFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                        const std::shared_ptr<Distmat>& distmat,
                        const int grain,
                        const int nrows)
        : ParallelWorker(grain, 10, 1000)
        , dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , nrows_(nrows)
    { }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        id_t i = nrows_;
        id_t j;
        for (id_t id = begin; id < end; id++) {
            if (is_interrupted(id)) break; // nocov

            if (i >= nrows_ - 1)
                s2d(id, nrows_, i, j); // move to next column
            else
                i++; // same column still valid, only increase row

            double dist = dist_calculator->calculate(i,j);
            this->update_distmat(id, i, j, dist);
        }

        // delete local calculator
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    const std::shared_ptr<DistanceCalculator> dist_calculator_;

protected:
    std::shared_ptr<Distmat> distmat_;

    virtual void update_distmat(id_t id, id_t i, id_t j, double dist) = 0;
    //     (*distmat_)(i,j) = dist;
    //     (*distmat_)(j,i) = dist;
    // }

private:
    const id_t nrows_;
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
// void SymmetricFiller::fill() const {
//     // number of elements in square matrix without including diagonal
//     int nrows = distmat_->nrow();
//     int size = nrows * (nrows - 1) / 2;
//     int grain = get_grain(size, num_threads_);
//     SymmetricFillWorker fill_worker(dist_calculator_, distmat_, grain, distmat_->nrow());
//     parallel_for(0, size, fill_worker, grain);
// }

// =================================================================================================
/* multi-threaded lower triangular filler
 *   This filler assumes a single column vector has been received in *distmat_
 */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LowerTriangularFiller::LowerTriangularFiller(std::shared_ptr<Distmat>& distmat,
                                             const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                             const SEXP& NUM_THREADS)
    : DistmatFiller(distmat, dist_calculator, NUM_THREADS)
{ }

// -------------------------------------------------------------------------------------------------
/* parallel worker */
// -------------------------------------------------------------------------------------------------
class LowerTriangularFillWorker : public SymmetricFillWorker {
public:
    // constructor
    LowerTriangularFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                              const std::shared_ptr<Distmat>& distmat,
                              const int grain,
                              const int nrows)
        : SymmetricFillWorker(dist_calculator, distmat, grain, nrows)
    { }

protected:
    void update_distmat(id_t id, id_t i, id_t j, double dist) override {
        (*distmat_)(id,0) = dist;
    }
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void LowerTriangularFiller::fill() const {
    int size = distmat_->nrow();
    if (size <= 0) return;

    int rows_for_worker = static_cast<int>((1 + std::sqrt(1 + 8.0 * size)) / 2);
    int grain = get_grain(size, num_threads_);
    LowerTriangularFillWorker fill_worker(dist_calculator_, distmat_, grain, rows_for_worker);
    parallel_for(0, size, fill_worker, grain);
}

// =================================================================================================
/* multi-threaded lower triangular and diagonal filler
 *   This filler assumes a single column vector has been received in *distmat_
 */
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/* constructor */
// -------------------------------------------------------------------------------------------------
LowerTriangularDiagonalFiller::LowerTriangularDiagonalFiller(std::shared_ptr<Distmat>& distmat,
                                                             const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                                             const SEXP& NUM_THREADS)
    : DistmatFiller(distmat, dist_calculator, NUM_THREADS)
{ }

// -------------------------------------------------------------------------------------------------
/* parallel worker */
// -------------------------------------------------------------------------------------------------
class LowerTriangularDiagonalFillWorker : public ParallelWorker {
public:
    // constructor
    LowerTriangularDiagonalFillWorker(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                      const std::shared_ptr<Distmat>& distmat,
                                      const int grain,
                                      const int nrows)
        : ParallelWorker(grain, 10, 1000)
        , dist_calculator_(dist_calculator)
        , distmat_(distmat)
        , nrows_(nrows)
    { }

    // parallel loop across specified range
    void work_it(id_t begin, id_t end) override {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();

        // fill distances
        id_t i, j;
        begin_indices(begin, nrows_, i, j);
        for (id_t id = begin; id < end; id++) {
            if (is_interrupted(id)) break; // nocov

            if (dist_calculator->distance == "SDTW" || i != j) {
                double dist = dist_calculator->calculate(i,j);
                (*distmat_)(id,0) = dist;
            }

            i++;
            if (i >= nrows_) {
                j++;
                i = j;
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
    const id_t nrows_;

    void begin_indices(const id_t id, const id_t nrow, id_t& i, id_t& j) {
        j = 0;
        id_t upper_bound = nrow - 1;
        while (id > upper_bound) {
            j++;
            upper_bound += nrow - j;
        }
        i = nrow - 1 - (upper_bound - id);
    }
};

// -------------------------------------------------------------------------------------------------
/* public fill method */
// -------------------------------------------------------------------------------------------------
void LowerTriangularDiagonalFiller::fill() const {
    int size = distmat_->nrow();
    if (size <= 0) return;

    int rows_for_worker = static_cast<int>((-1 + std::sqrt(1 + 8.0 * size)) / 2);
    int grain = get_grain(size, num_threads_);
    LowerTriangularDiagonalFillWorker fill_worker(dist_calculator_, distmat_, grain, rows_for_worker);
    parallel_for(0, size, fill_worker, grain);
}

} // namespace dtwclust
