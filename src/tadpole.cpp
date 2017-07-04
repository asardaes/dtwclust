#include <Rcpp.h>
#include <algorithm> // std::stable_sort and std::sort
#include <iomanip> // std::setprecision
#include <vector>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* class that stores lower triangular of a matrix and knows how to access it */
// =================================================================================================

template <typename T>
class LowerTriMat {
public:
    // ---------------------------------------------------------------------------------------------
    /* constructors */
    // ---------------------------------------------------------------------------------------------

    LowerTriMat(int size) : _size(size)
    {
        if (size < 1) Rcpp::stop("TADPole: invalid dimension for a distance matrix");
        _len = size * (size + 1) / 2 - size;
        _data = new T[_len];
        for (int i = 0; i < _len; i++) _data[i] = T(0);
    }

    LowerTriMat(int size, T init_val) : _size(size)
    {
        if (size < 1) Rcpp::stop("TADPole: invalid dimension for a distance matrix");
        _len = size * (size + 1) / 2 - size;
        _data = new T[_len];
        for (int i = 0; i < _len; i++) _data[i] = init_val;
    }

    // ---------------------------------------------------------------------------------------------
    /* copy, assign and destructor (rule of 3) */
    // ---------------------------------------------------------------------------------------------

    LowerTriMat(const LowerTriMat& ltm)
    {
        _size = ltm._size;
        _len = ltm._len;
        _data = new T[_len];
        for (int i = 0; i < _len; i++) _data[i] = ltm._data[i];
    }

    LowerTriMat& operator= (const LowerTriMat& ltm)
    {
        _size = ltm._size;
        _len = ltm._len;
        _data = new T[_len];
        for (int i = 0; i < _len; i++) _data[i] = ltm._data[i];
        return *this;
    }

    ~LowerTriMat() { delete[] _data; }

    // ---------------------------------------------------------------------------------------------
    /* operator() */
    // ---------------------------------------------------------------------------------------------

    T& operator() (int row, int col)
    {
        if (row >= _size || col >= _size || row == col)
            Rcpp::stop("TADPole: invalid indices for a distance matrix");
        if (col > row) {
            int swap = row;
            row = col;
            col = swap;
        }
        return _data[d2s(row, col)];
    }

    const T operator() (int row, int col) const
    {
        if (row >= _size || col >= _size || row == col)
            Rcpp::stop("TADPole: invalid indices for a distance matrix");
        if (col > row) {
            int swap = row;
            row = col;
            col = swap;
        }
        return _data[d2s(row, col)];
    }

private:
    int _size, _len;
    T* _data;

    // double to single index with adjustment for missing upper triangular
    int d2s(const int row, const int col) const
    {
        int adjustment = 0;
        for (int k = col; k >= 0; k--) adjustment += k + 1;
        int id = row + col*_size - adjustment;
        if (id >= _len)
            Rcpp::stop("Something went wrong, an invalid distance matrix index was computed");
        return id;
    }
};

// =================================================================================================
/* pruning during local density calculation */
// =================================================================================================

std::vector<double> local_density(const Rcpp::List& series,
                                  const int num_series,
                                  double dc,
                                  const Rcpp::List& dtw_args,
                                  const Rcpp::NumericMatrix& LBM,
                                  const Rcpp::NumericMatrix& UBM,
                                  LowerTriMat<double>& distmat,
                                  LowerTriMat<int>& flags,
                                  int& num_dist_op)
{
    std::vector<double> rho(num_series, 0);

    /*
     * Flag definition
     *   0 - DTW calculated, and it lies below dc
     *   1 - calculate DTW
     *   2 - within dc, prune
     *   3 - not within dc, prune
     *   4 - identical series
     */
    for (int i = 1; i < num_series; i++) {
        R_CheckUserInterrupt();
        for (int j = 0; j < i; j++) {
            if (LBM(i,j) <= dc && UBM(i,j) > dc) {
                num_dist_op++;
                Rcpp::NumericVector x = series[i];
                Rcpp::NumericVector y = series[j];
                double dtw_dist = dtwb(x, y, dtw_args);
                distmat(i,j) = dtw_dist;
                if (dtw_dist <= dc)
                    flags(i,j) = 0;
                else
                    flags(i,j) = 1;

            } else if (LBM(i,j) <= dc && UBM(i,j) < dc) {
                flags(i,j) = 2;
            } else if (LBM(i,j) > dc) {
                flags(i,j) = 3;
            } else {
                flags(i,j) = 4;
            }
        }
    }

    bool no_peaks = true;
    for (int i = 0; i < num_series; i++) {
        for (int j = 0; j < num_series; j++) {
            if (i == j) continue;
            int flag = flags(i,j);
            if (flag == 0 || flag == 2) rho[i] += 1;
        }
        if (rho[i] > 0) no_peaks = false;
    }

    if (no_peaks)
        Rcpp::stop("No density peaks detected, choose a different value for cutoff distance 'dc'");

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
            this_rho = 1;
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
        R_CheckUserInterrupt();
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
    reorder(delta_ub, id_orig); // template in dtwclustpp.h
    return delta_ub;
}

// =================================================================================================
/* pruning during NN distance calculation from higher density list (phase 2) */
// =================================================================================================

std::vector<double> nn_dist_2(const Rcpp::List& series,
                              const int num_series,
                              const Rcpp::List& dtw_args,
                              const std::vector<size_t>& id_cl,
                              const std::vector<double>& delta_ub,
                              const Rcpp::NumericMatrix& LBM,
                              const Rcpp::NumericMatrix& UBM,
                              const LowerTriMat<int>& flags,
                              const LowerTriMat<double>& distmat,
                              std::vector<int>& nearest_neighbors,
                              int& num_dist_op)
{
    std::vector<double> delta(num_series);
    nearest_neighbors[0] = -1;

    double min_delta = R_PosInf, max_delta = 0;
    for (int i = 1; i < num_series; i++) {
        int which_min_delta = -1;
        min_delta = R_PosInf;
        for (int j = 0; j < i; j++) {
            R_CheckUserInterrupt();
            int ii = id_cl[i], jj = id_cl[j];
            bool prune = LBM(ii,jj) > delta_ub[ii];
            bool precomputed = flags(ii,jj) == 0 || flags(ii,jj) == 1;
            if (precomputed) {
                double dtw_dist = distmat(ii,jj);
                if (dtw_dist < min_delta) {
                    min_delta = dtw_dist;
                    which_min_delta = jj;
                }
            } else if (prune) {
                double ub = UBM(ii,jj);
                if (ub < min_delta) {
                    min_delta = ub;
                    which_min_delta = jj;
                }
            } else {
                num_dist_op++;
                Rcpp::NumericVector x = series[ii];
                Rcpp::NumericVector y = series[jj];
                double dtw_dist = dtwb(x, y, dtw_args);
                if (dtw_dist < min_delta) {
                    min_delta = dtw_dist;
                    which_min_delta = jj;
                }
            }
        }
        if (min_delta > max_delta) max_delta = min_delta;
        delta[i] = min_delta;
        nearest_neighbors[i] = which_min_delta;
    }
    delta[0] = max_delta;

    // get min and max
    min_delta = num_series + 1;
    max_delta = -1;
    for (double const &this_delta : delta) {
        if (this_delta < min_delta) min_delta = this_delta;
        if (this_delta > max_delta) max_delta = this_delta;
    }

    // normalize
    double den = max_delta - min_delta;
    for (double& this_delta : delta) {
        if (den == 0)
            this_delta = 1;
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
                        std::vector<size_t>& id_cent,
                        const std::vector<size_t>& id_cl,
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
            int ii = (int)(id_cent[i]);
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
            Rcpp::Rcout << "TADPole completed for k = " << k << " & dc = " << dc << "\n\n";

        list(counter) = Rcpp::List::create(
            Rcpp::_["cl"] = cl,
            Rcpp::_["centroids"] = cent,
            Rcpp::_["distCalcPercentage"] = dist_op_percent
        );
    }
}

// =================================================================================================
/* main C++ function */
// =================================================================================================

SEXP tadpole_cpp(const Rcpp::List& series,
                 const Rcpp::IntegerVector& k,
                 const double dc,
                 const Rcpp::List& dtw_args,
                 const Rcpp::NumericMatrix& LBM,
                 const Rcpp::NumericMatrix& UBM,
                 const bool trace,
                 Rcpp::List& list)
{
    int num_series = series.length();
    LowerTriMat<double> distmat(num_series, NA_REAL);
    LowerTriMat<int> flags(num_series, -1);
    int num_dist_op = 0;

    if (trace) Rcpp::Rcout << "\tPruning during local density calculation\n";
    Rflush();
    std::vector<double> rho = local_density(series, num_series,
                                            dc, dtw_args,
                                            LBM, UBM,
                                            distmat, flags, num_dist_op);

    if (trace) Rcpp::Rcout << "\tPruning during nearest-neighbor distance calculation (phase 1)\n";
    Rflush();
    std::vector<double> delta_ub = nn_dist_1(rho, num_series, distmat, UBM);

    // get indices of sorted rho (using template function from dtwclustpp.h)
    std::vector<double> helper = rho;
    for (int i = 0; i < num_series; i++) helper[i] *= delta_ub[i];
    auto id_cl = stable_sort_ind(helper, true);

    if (trace) Rcpp::Rcout << "\tPruning during nearest-neighbor distance calculation (phase 2)\n";
    Rflush();
    std::vector<int> nearest_neighbors(num_series);
    std::vector<double> delta = nn_dist_2(series, num_series,
                                          dtw_args, id_cl, delta_ub,
                                          LBM, UBM, flags, distmat,
                                          nearest_neighbors, num_dist_op);

    auto id_orig = stable_sort_ind(id_cl, false);
    reorder(delta, id_orig);
    for (int i = 0; i < num_series; i++) helper[i] = rho[i] * delta[i];
    auto id_cent = stable_sort_ind(helper, true);
    double dist_op_percent = (num_dist_op / ((double)num_series * (num_series + 1) / 2 - num_series)) * 100;

    if (trace) {
        Rcpp::Rcout << "\tPruning percentage = " << std::setprecision(3) << 100 - dist_op_percent << "%\n";
        Rcpp::Rcout << "\tPerforming cluster assignmnet\n\n";
        Rflush();
    }
    cluster_assignment(k, dc, id_cent, id_cl, nearest_neighbors, dist_op_percent, trace, list);

    return R_NilValue;
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
                        SEXP LB, SEXP UB, SEXP TRACE,
                        SEXP LIST)
{
    BEGIN_RCPP
    Rcpp::List series(X), dtw_args(DTW_ARGS), list(LIST);
    Rcpp::NumericMatrix LBM(LB), UBM(UB);
    Rcpp::IntegerVector k(K);
    double dc = Rcpp::as<double>(DC);
    bool trace = Rcpp::as<bool>(TRACE);
    return tadpole_cpp(series, k, dc, dtw_args, LBM, UBM, trace, list);
    END_RCPP
}

} // namespace dtwclust
