#include "ParallelWorker.h"

#include <cstddef> // std::size_t

#include <RcppParallel.h>
#undef ERROR // collision between R.h and mingw_32/i686-w64-mingw32/include/windows.h
#include <RcppThread.h>

namespace dtwclust {

void ParallelWorker::operator()(std::size_t begin, std::size_t end) {
    work_it(begin, end);
    // make sure this is called at least once per thread call
    RcppThread::isInterrupted();
}

ParallelWorker::ParallelWorker(const int grain, const int min, const int max)
    : interrupt_grain_(interrupt_grain(grain, min, max))
{ }

bool ParallelWorker::is_interrupted() const {
    return RcppThread::isInterrupted();
}

bool ParallelWorker::is_interrupted(const std::size_t i) const {
    return RcppThread::isInterrupted(i % interrupt_grain_ == 0);
}

// how often to check for user interrupt inside a thread
int ParallelWorker::interrupt_grain(const int grain, const int min, const int max) const {
    int result = grain / 1000;
    if (result < min) result = min;
    if (result > max) result = max;
    if (result < 1) result = 1;
    return result;
}

} // namespace dtwclust
