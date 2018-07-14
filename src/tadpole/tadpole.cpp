#include "R-gateways.h"

#include <algorithm> // std::sort
#include <atomic> // atomic_int
#include <cstddef> // std::size_t
#include <iomanip> // std::setprecision
#include <memory> // *_ptr
#include <string>
#include <vector>

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "../distances/calculators.h"
#include "../utils/utils.h" // Rflush, get_grain, s2d, id_t

namespace dtwclust {

// =================================================================================================
/* Templates */
// =================================================================================================

// see https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> stable_sort_ind(const std::vector<T>& v, const bool decreasing)
{
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    if (decreasing)
        std::stable_sort(idx.begin(), idx.end(),
                         [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
    else
        std::stable_sort(idx.begin(), idx.end(),
                         [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}

// see https://stackoverflow.com/a/22183350/5793905
template <typename T>
void reorder(std::vector<T>& v, std::vector<size_t>& order)
{
    // for all elements to put in place
    for (size_t i = 0; i < v.size(); ++i) {
        // while order[i] is not yet in place
        // every swap places at least one element in its proper place
        while (order[i] != order[order[i]]) {
            std::swap(v[order[i]], v[order[order[i]]]);
            std::swap(order[i], order[order[i]]);
        }
    }
}

// =================================================================================================
/* class that stores lower triangular of a matrix and knows how to access it */
// =================================================================================================

template <typename T>
class LowerTriMat {
public:
    // ---------------------------------------------------------------------------------------------
    /* constructors */
    // ---------------------------------------------------------------------------------------------

    LowerTriMat(int nrow) : nrow_(nrow)
    {
        if (nrow < 1) Rcpp::stop("TADPole: invalid dimension for a distance matrix"); // nocov
        len_ = nrow * (nrow + 1) / 2 - nrow;
        data_ = new T[len_];
        for (int i = 0; i < len_; i++) data_[i] = T(0);
    }

    LowerTriMat(int nrow, T init_val) : nrow_(nrow)
    {
        if (nrow < 1) Rcpp::stop("TADPole: invalid dimension for a distance matrix"); // nocov
        len_ = nrow * (nrow + 1) / 2 - nrow;
        data_ = new T[len_];
        for (int i = 0; i < len_; i++) data_[i] = init_val;
    }

    // ---------------------------------------------------------------------------------------------
    /* copy, assign and destructor (rule of 3) */
    // ---------------------------------------------------------------------------------------------

    LowerTriMat(const LowerTriMat& ltm)
    {
        nrow_ = ltm.nrow_;
        len_ = ltm.len_;
        data_ = new T[len_];
        for (int i = 0; i < len_; i++) data_[i] = ltm.data_[i];
    }

    LowerTriMat& operator= (const LowerTriMat& ltm)
    {
        nrow_ = ltm.nrow_;
        len_ = ltm.len_;
        data_ = new T[len_];
        for (int i = 0; i < len_; i++) data_[i] = ltm.data_[i];
        return *this;
    }

    ~LowerTriMat() { delete[] data_; }

    // ---------------------------------------------------------------------------------------------
    /* operator() */
    // ---------------------------------------------------------------------------------------------

    T& operator()(int row, int col)
    {
        if (row >= nrow_ || col >= nrow_ || row == col)
            Rcpp::stop("TADPole: invalid indices for a distance matrix"); // nocov
        if (col > row) {
            int swap = row;
            row = col;
            col = swap;
        }
        return data_[d2s(row, col)];
    }

    const T operator()(int row, int col) const
    {
        if (row >= nrow_ || col >= nrow_ || row == col)
            Rcpp::stop("TADPole: invalid indices for a distance matrix"); // nocov
        if (col > row) {
            int swap = row;
            row = col;
            col = swap;
        }
        return data_[d2s(row, col)];
    }

    // ---------------------------------------------------------------------------------------------
    /* operator[] */
    // ---------------------------------------------------------------------------------------------

    T& operator[](const int id) { return data_[id]; }
    const T operator[](const int id) const { return data_[id]; }

    // ---------------------------------------------------------------------------------------------
    /* others */
    // ---------------------------------------------------------------------------------------------

    int length() { return len_; }

    // ---------------------------------------------------------------------------------------------
    /* private members */
    // ---------------------------------------------------------------------------------------------
private:
    int nrow_, len_;
    T* data_;

    // double to single index with adjustment for missing upper triangular and diagonal
    int d2s(const int row, const int col) const {
        int adjustment = 0;
        for (int k = col; k >= 0; k--) adjustment += k + 1;
        int id = row + col*nrow_ - adjustment;
        return id;
    }
};

// =================================================================================================
/* helper classes for multi-threading */
// =================================================================================================

// local density
class LocalDensityHelper : public RcppParallel::Worker {
public:
    // constructor
    LocalDensityHelper(const double dc,
                       const std::shared_ptr<DistanceCalculator>& dist_calculator,
                       const Rcpp::NumericMatrix& LBM,
                       const Rcpp::NumericMatrix& UBM,
                       LowerTriMat<double>& distmat,
                       LowerTriMat<int>& flags,
                       std::atomic_int& num_dist_op)
        : dc_(dc)
        , dist_calculator_(dist_calculator)
        , LBM_(LBM)
        , UBM_(UBM)
        , distmat_(distmat)
        , flags_(flags)
        , num_dist_op_(num_dist_op)
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();
        /*
         * Flag definition
         *   0 - DTW calculated, and it lies below dc
         *   1 - calculated DTW, but it's above dc
         *   2 - within dc, prune
         *   3 - not within dc, prune
         *   4 - identical series
         */
        id_t i, j;
        for (std::size_t id = begin; id < end; id++) {
            s2d(id, LBM_.nrow(), i, j);
            if (LBM_(i,j) <= dc_ && UBM_(i,j) > dc_) {
                num_dist_op_++;
                double dtw_dist = dist_calculator->calculate(i,j);
                distmat_[id] = dtw_dist;
                if (dtw_dist <= dc_)
                    flags_[id] = 0;
                else
                    flags_[id] = 1;
            }
            else if (LBM_(i,j) <= dc_ && UBM_(i,j) < dc_) {
                flags_[id] = 2;
            }
            else if (LBM_(i,j) > dc_) {
                flags_[id] = 3;
            }
            else {
                flags_[id] = 4; // nocov
            }
        }
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    const double dc_;
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    const RcppParallel::RMatrix<double> LBM_, UBM_;
    LowerTriMat<double>& distmat_;
    LowerTriMat<int>& flags_;
    std::atomic_int& num_dist_op_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// pruning during NN calculation (phase 2)
class PruningHelper : public RcppParallel::Worker {
public:
    // constructor
    PruningHelper(const std::shared_ptr<DistanceCalculator>& dist_calculator,
                  const std::vector<std::size_t>& id_cl,
                  const std::vector<double>& delta_ub,
                  const Rcpp::NumericMatrix& LBM,
                  const Rcpp::NumericMatrix& UBM,
                  const LowerTriMat<int>& flags,
                  const LowerTriMat<double>& distmat,
                  std::vector<int>& nearest_neighbors,
                  std::vector<double>& delta,
                  std::atomic_int& num_dist_op,
                  double& max_delta)
        : dist_calculator_(dist_calculator)
        , LBM_(LBM)
        , UBM_(UBM)
        , flags_(flags)
        , distmat_(distmat)
        , id_cl_(id_cl)
        , delta_ub_(delta_ub)
        , delta_(delta)
        , nearest_neighbors_(nearest_neighbors)
        , num_dist_op_(num_dist_op)
        , max_delta_(max_delta)
    { }

    // parallel loop across specified range
    void operator()(std::size_t begin, std::size_t end) {
        // local copy of dist_calculator so it is setup separately for each thread
        mutex_.lock();
        DistanceCalculator* dist_calculator = dist_calculator_->clone();
        mutex_.unlock();
        for (std::size_t i = begin; i < end; i++) {
            int which_min_delta = -1;
            double min_delta = R_PosInf;
            for (std::size_t j = 0; j < i; j++) {
                std::size_t ii = id_cl_[i], jj = id_cl_[j];
                bool prune = LBM_(ii,jj) > delta_ub_[ii];
                bool precomputed = flags_(ii,jj) == 0 || flags_(ii,jj) == 1;
                if (precomputed) {
                    double dtw_dist = distmat_(ii,jj);
                    if (dtw_dist < min_delta) {
                        min_delta = dtw_dist;
                        which_min_delta = jj;
                    }
                }
                else if (prune) {
                    double ub = UBM_(ii,jj);
                    if (ub < min_delta) {
                        min_delta = ub;
                        which_min_delta = jj;
                    }
                }
                else {
                    num_dist_op_++;
                    double dtw_dist = dist_calculator->calculate(ii,jj);
                    if (dtw_dist < min_delta) {
                        min_delta = dtw_dist;
                        which_min_delta = jj;
                    }
                }
            }
            delta_[i] = min_delta;
            nearest_neighbors_[i] = which_min_delta;
            mutex_.lock();
            if (min_delta > max_delta_) max_delta_ = min_delta;
            mutex_.unlock();
        }
        mutex_.lock();
        delete dist_calculator;
        mutex_.unlock();
    }

private:
    const std::shared_ptr<DistanceCalculator> dist_calculator_;
    const RcppParallel::RMatrix<double> LBM_, UBM_;
    const LowerTriMat<int>& flags_;
    const LowerTriMat<double>& distmat_;
    const std::vector<std::size_t>& id_cl_;
    const std::vector<double>& delta_ub_;
    std::vector<double>& delta_;
    std::vector<int>& nearest_neighbors_;
    std::atomic_int& num_dist_op_;
    double& max_delta_;
    // for synchronization during memory allocation (from TinyThread++, comes with RcppParallel)
    tthread::mutex mutex_;
};

// =================================================================================================
/* pruning during local density calculation */
// =================================================================================================

std::vector<double> local_density(const Rcpp::List& series,
                                  const int num_series,
                                  const double dc,
                                  const std::shared_ptr<DistanceCalculator>& dist_calculator,
                                  const Rcpp::NumericMatrix& LBM,
                                  const Rcpp::NumericMatrix& UBM,
                                  LowerTriMat<double>& distmat,
                                  LowerTriMat<int>& flags,
                                  std::atomic_int& num_dist_op,
                                  const int num_threads)
{
    std::vector<double> rho(num_series, 0);
    int grain = get_grain(distmat.length(), num_threads);
    LocalDensityHelper parallel_worker(
            dc,
            dist_calculator,
            LBM,
            UBM,
            distmat,
            flags,
            num_dist_op
    );
    RcppParallel::parallelFor(0, distmat.length(), parallel_worker, grain);
    bool no_peaks = true;
    for (int i = 0; i < num_series; i++) {
        for (int j = 0; j < num_series; j++) {
            if (i == j) continue;
            int flag = flags(i,j);
            if (flag == 0 || flag == 2) rho[i] += 1;
        }
        if (rho[i] > 0) no_peaks = false;
    }
    if (no_peaks) Rcpp::stop("No density peaks detected, choose a different value for cutoff distance 'dc'");
    // get min and max
    double min_rho = num_series+1, max_rho = -1;
    for (double const &this_rho : rho) {
        if (this_rho < min_rho) min_rho = this_rho;
        if (this_rho > max_rho) max_rho = this_rho;
    }
    // normalize
    double den = max_rho - min_rho;
    for (double& this_rho : rho) {
        if (den == 0)
            this_rho = 1; // nocov
        else
            this_rho = (this_rho - min_rho) / den;
    }
    return rho;
}

// =================================================================================================
/* pruning during NN distance calculation from higher density list (phase 1) */
// =================================================================================================

std::vector<double> nn_dist_1(const std::vector<double>& rho, const int num_series,
                              const LowerTriMat<double>& distmat, const Rcpp::NumericMatrix& UBM)
{
    std::vector<double> delta_ub(num_series);
    auto id_rho_sorted = stable_sort_ind(rho, true);
    double max_delta = 0;
    for (int i = 1; i < num_series; i++) {
        double min_ub_i = R_PosInf;
        for (int j = 0; j < i; j++) {
            int ii = id_rho_sorted[i], jj = id_rho_sorted[j];
            double ub_i = distmat(ii,jj);
            if (Rcpp::NumericVector::is_na(ub_i)) ub_i = UBM(ii,jj);
            if (ub_i < min_ub_i) min_ub_i = ub_i;
        }
        delta_ub[i] = min_ub_i;
        if (min_ub_i > max_delta) max_delta = min_ub_i;
    }
    delta_ub[0] = max_delta;
    auto id_orig = stable_sort_ind(id_rho_sorted, false);
    reorder(delta_ub, id_orig); // template
    return delta_ub;
}

// =================================================================================================
/* pruning during NN distance calculation from higher density list (phase 2) */
// =================================================================================================

std::vector<double> nn_dist_2(const Rcpp::List& series,
                              const int num_series,
                              const std::shared_ptr<DistanceCalculator>& dist_calculator,
                              const std::vector<std::size_t>& id_cl,
                              const std::vector<double>& delta_ub,
                              const Rcpp::NumericMatrix& LBM,
                              const Rcpp::NumericMatrix& UBM,
                              const LowerTriMat<int>& flags,
                              const LowerTriMat<double>& distmat,
                              std::vector<int>& nearest_neighbors,
                              std::atomic_int& num_dist_op,
                              const int num_threads)
{
    std::vector<double> delta(num_series);
    nearest_neighbors[0] = -1;
    double max_delta = 0;
    int grain = get_grain(num_series, num_threads);
    PruningHelper parallel_worker(
            dist_calculator,
            id_cl,
            delta_ub,
            LBM,
            UBM,
            flags,
            distmat,
            nearest_neighbors,
            delta,
            num_dist_op,
            max_delta
    );
    RcppParallel::parallelFor(1, num_series, parallel_worker, grain);
    delta[0] = max_delta;
    // get min and max
    double min_delta = num_series + 1;
    max_delta = -1;
    for (const double& this_delta : delta) {
        if (this_delta < min_delta) min_delta = this_delta;
        if (this_delta > max_delta) max_delta = this_delta;
    }
    // normalize
    double den = max_delta - min_delta;
    for (double& this_delta : delta) {
        if (den == 0)
            this_delta = 1; // nocov
        else
            this_delta = (this_delta - min_delta) / den;
    }
    return delta;
}

// =================================================================================================
/* cluster assignment */
// =================================================================================================

void cluster_assignment(const Rcpp::IntegerVector& k_vec,
                        const double dc,
                        std::vector<std::size_t>& id_cent,
                        const std::vector<std::size_t>& id_cl,
                        const std::vector<int>& nearest_neighbors,
                        const double dist_op_percent,
                        const bool trace,
                        Rcpp::List& list)
{
    int len = k_vec.length();
    for (int counter = 0; counter < len; counter++) {
        int k = k_vec[counter];
        int n = id_cl.size();
        Rcpp::IntegerVector cl = Rcpp::rep(-1, n); // cluster ids
        Rcpp::IntegerVector cent(k); // centroid ids
        // id_cent only contains distinct elements, so stable sorting is not needed
        std::sort(id_cent.begin(), id_cent.begin() + k);
        for (int i = 0; i < k; i++) {
            int ii = static_cast<int>(id_cent[i]);
            cent[i] = ii + 1;
            cl[ii] = i + 1;
        }
        bool warn = false;
        for (int i = 0; i < n; i++) {
            int ii = id_cl[i];
            if (cl[ii] == -1) {
                cl[ii] = cl[nearest_neighbors[i]];
                if (cl[ii] == -1) warn = true;
            }
        }
        if (warn) // nocov start
            Rcpp::warning(
                "At least one series wasn't assigned to a cluster. This shouldn't happen, please contact maintainer."
            ); // nocov end
        if (trace)
            Rcpp::Rcout << "TADPole completed for k = " << k << " & dc = " << dc << "\n";
        list(counter) = Rcpp::List::create(
            Rcpp::_["cl"] = cl,
            Rcpp::_["centroids"] = cent,
            Rcpp::_["distCalcPercentage"] = dist_op_percent
        );
    }
    if (trace) Rcpp::Rcout << "\n";
}

// =================================================================================================
/* main C++ function */
// =================================================================================================

SEXP tadpole_cpp(const Rcpp::List& series,
                 const Rcpp::IntegerVector& k,
                 const double dc,
                 const SEXP& DTW_ARGS,
                 const Rcpp::NumericMatrix& LBM,
                 const Rcpp::NumericMatrix& UBM,
                 const bool trace,
                 Rcpp::List& list,
                 const int num_threads)
{
    auto dist_calculator = DistanceCalculatorFactory().create(
        "DTW_BASIC", DTW_ARGS, series, series);

    int num_series = series.length();
    LowerTriMat<double> distmat(num_series, NA_REAL);
    LowerTriMat<int> flags(num_series, -1);
    std::atomic_int num_dist_op(0);

    if (trace) Rcpp::Rcout << "Pruning during local density calculation\n";
    Rflush();
    std::vector<double> rho = local_density(
        series, num_series, dc, dist_calculator,
        LBM, UBM, distmat, flags, num_dist_op, num_threads);

    if (trace) Rcpp::Rcout << "Pruning during nearest-neighbor distance calculation (phase 1)\n";
    Rflush();
    std::vector<double> delta_ub = nn_dist_1(rho, num_series, distmat, UBM);

    // get indices of sorted rho
    std::vector<double> helper = rho;
    for (int i = 0; i < num_series; i++) helper[i] *= delta_ub[i];
    auto id_cl = stable_sort_ind(helper, true);

    if (trace) Rcpp::Rcout << "Pruning during nearest-neighbor distance calculation (phase 2)\n";
    Rflush();
    std::vector<int> nearest_neighbors(num_series);
    std::vector<double> delta = nn_dist_2(
        series, num_series, dist_calculator, id_cl, delta_ub,
        LBM, UBM, flags, distmat, nearest_neighbors, num_dist_op, num_threads);

    auto id_orig = stable_sort_ind(id_cl, false);
    reorder(delta, id_orig);
    for (int i = 0; i < num_series; i++) helper[i] = rho[i] * delta[i];
    auto id_cent = stable_sort_ind(helper, true);
    double dist_op_percent = (num_dist_op / ((double)num_series * (num_series + 1) / 2 - num_series)) * 100;

    if (trace) {
        Rcpp::Rcout << "Pruning percentage = " <<
            std::setprecision(3) << 100 - dist_op_percent << "%\n";
        Rcpp::Rcout << "Performing cluster assignment\n\n";
        Rflush();
    }
    cluster_assignment(k, dc, id_cent, id_cl, nearest_neighbors, dist_op_percent, trace, list);
    return R_NilValue;
}

// =================================================================================================
/* gateway function */
// =================================================================================================

extern "C" SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST, SEXP NUM_THREADS)
{
    BEGIN_RCPP
    Rcpp::List list(LIST);
    Rcpp::NumericMatrix LBM(LB), UBM(UB);
    Rcpp::IntegerVector k(K);
    double dc = Rcpp::as<double>(DC);
    bool trace = Rcpp::as<bool>(TRACE);
    int num_threads = Rcpp::as<int>(NUM_THREADS);
    return tadpole_cpp(X, k, dc, DTW_ARGS, LBM, UBM, trace, list, num_threads);
    END_RCPP
}

} // namespace dtwclust
