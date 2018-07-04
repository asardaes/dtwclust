#include "R-gateways.h"

#include <cmath> // std::round
#include <unordered_set>
#include <vector>

#include <Rcpp.h>

#include "UndirectedGraph.h"
#include "utils.h"

namespace dtwclust {

// =================================================================================================
/* Helper for SparseDistmat in R */
// =================================================================================================

class SparseDistmatIndices
{
public:
    SparseDistmatIndices(int num_rows)
        : num_rows_(num_rows)
    { }

    Rcpp::IntegerMatrix getNewIndices(
            const Rcpp::IntegerVector& i,
            const Rcpp::IntegerVector& j,
            const bool symmetric)
    {
        std::vector<int> new_i, new_j;
        for (int ii = 0; ii < i.length(); ii++) {
            for (int jj = 0; jj < j.length(); jj++) {
                int this_i = i[ii];
                int this_j = j[jj];
                if (symmetric && this_j > this_i) {
                    int swap = this_i;
                    this_i = this_j;
                    this_j = swap;
                }
                int ij = this_i + (this_j - 1) * num_rows_;
                bool new_index = existing_indices_.find(ij) == existing_indices_.end();
                if (new_index) {
                    existing_indices_.insert(ij);
                    new_i.push_back(this_i);
                    new_j.push_back(this_j);
                }
            }
        }
        Rcpp::IntegerMatrix new_ids(new_i.size(), 2);
        for (int k = 0; k < new_ids.nrow(); k++) {
            new_ids(k, 0) = new_i[k];
            new_ids(k, 1) = new_j[k];
        }
        return new_ids;
    }

private:
    int num_rows_;
    std::unordered_set<int> existing_indices_;
};

extern "C" SEXP SparseDistmatIndices__new(SEXP num_rows)
{
    BEGIN_RCPP
    Rcpp::XPtr<SparseDistmatIndices> ptr(new SparseDistmatIndices(Rcpp::as<int>(num_rows)), true);
    return ptr;
    END_RCPP
}

extern "C" SEXP SparseDistmatIndices__getNewIndices(SEXP xptr, SEXP i, SEXP j, SEXP symmetric)
{
    BEGIN_RCPP
    Rcpp::XPtr<SparseDistmatIndices> ptr(xptr);
    return ptr->getNewIndices(i, j, Rcpp::as<bool>(symmetric));
    END_RCPP
}

// =================================================================================================
/* Gateway for warping envelopes */
// ================================================================================================

extern "C" SEXP envelope(SEXP series, SEXP window) {
    BEGIN_RCPP
    Rcpp::NumericVector x(series);
    int n = x.length();
    Rcpp::NumericVector L(n), U(n);
    SurrogateMatrix<double> temp_x(n, 1, &x[0]);
    SurrogateMatrix<double> temp_l(n, 1, &L[0]);
    SurrogateMatrix<double> temp_u(n, 1, &U[0]);
    envelope_cpp(temp_x, Rcpp::as<unsigned int>(window), temp_l, temp_u);
    Rcpp::List ret;
    ret["lower"] = L;
    ret["upper"] = U;
    return(ret);
    END_RCPP
}

// =================================================================================================
/* Force symmetry helper */
// =================================================================================================

extern "C" SEXP force_lb_symmetry(SEXP X)
{
    BEGIN_RCPP
    Rcpp::NumericMatrix matrix(X);
    for (int i = 1; i < matrix.nrow(); i++) {
        Rcpp::checkUserInterrupt();
        for (int j = 0; j < i; j++) {
            double lb1 = matrix(i,j);
            double lb2 = matrix(j,i);
            if (lb1 > lb2)
                matrix(j,i) = lb1;
            else
                matrix(i,j) = lb2;
        }
    }
    return R_NilValue;
    END_RCPP
}

// =================================================================================================
/* all possible combinations in pairs */
// =================================================================================================

void pairs_c(const int n, const int nrow, int *out)
{
    int i, j;
    int p = 0;
    for(j = 1; j < n; j++)
    {
        for(i = j + 1; i <= n; i++)
        {
            out[p] = i;
            out[p+nrow] = j;
            p++;
        }
    }
}

// the gateway function
extern "C" SEXP pairs(SEXP L)
{
    int n = Rf_asInteger(L);
    int nrow = n * (n+1) / 2 - n;

    // allocate output integer vector
    SEXP ret = PROTECT(Rf_allocMatrix(INTSXP, nrow, 2));

    // dispatch to C function
    pairs_c(n, nrow, INTEGER(ret));

    // release protection
    UNPROTECT(1);

    // finish
    return ret;
}

// =================================================================================================
/* assign existing names to existing vector */
// =================================================================================================

extern "C" SEXP setnames_inplace(SEXP vec, SEXP names) {
    Rf_setAttrib(vec, R_NamesSymbol, names);
    return R_NilValue;
}

// =================================================================================================
/* Helper for PairTracker in R */
// =================================================================================================

#define DONT_KNOW -1
#define CANNOT_LINK 0
#define MUST_LINK 1

/** Pair tracker
 *
 *  I could depend only on the aggregate_ graph without caring about the link_type, but then I
 *  wouldn't be able to give fine-grained feedback on the R side.
 */
class PairTracker
{
public:
    PairTracker(const int max_size)
        : must_link_(static_cast<unsigned int>(max_size))
        , cannot_link_(static_cast<unsigned int>(max_size))
        , dont_know_(static_cast<unsigned int>(max_size))
        , aggregate_(static_cast<unsigned int>(max_size))
        , max_size_(max_size)
    { }

    SEXP link(const int i, const int j, const int link_type) {
        if (i < 1 || i > max_size_ || j < 1 || j > max_size_) Rcpp::stop("Invalid indices provided");
        switch(link_type)
        {
        case DONT_KNOW:
            aggregate_.linkVertices(i,j);
            dont_know_.linkVertices(i,j);
            return Rcpp::wrap(dont_know_.isComplete());
        case CANNOT_LINK:
            aggregate_.linkVertices(i,j);
            cannot_link_.linkVertices(i,j);
            return Rcpp::wrap(cannot_link_.isComplete());
        case MUST_LINK:
            aggregate_.linkVertices(i,j, true);
            must_link_.linkVertices(i,j);
            return Rcpp::wrap(must_link_.isConnected());
        }
        return R_NilValue; // nocov
    }

    SEXP getUnseenPair() {
        if (aggregate_.isComplete()) return R_NilValue;
        Rcpp::IntegerVector pair(2);
        GetRNGstate(); // see https://github.com/rstudio/shiny/issues/1953
        bool seen = true;
        while (seen) {
            Rcpp::checkUserInterrupt();
            pair[0] = std::round(R::runif(1, max_size_));
            pair[1] = std::round(R::runif(1, max_size_));
            while (pair[0] == pair[1]) pair[1] = std::round(R::runif(1, max_size_));
            if (!aggregate_.areNeighbors(pair[0], pair[1]))
                seen = false;
        }
        PutRNGstate();
        return pair;
    }

private:
    UndirectedGraph must_link_, cannot_link_, dont_know_, aggregate_;
    int max_size_;
};

extern "C" SEXP PairTracker__new(SEXP max_size)
{
    BEGIN_RCPP
    Rcpp::XPtr<PairTracker> ptr(new PairTracker(Rcpp::as<int>(max_size)), true);
    return ptr;
    END_RCPP
}

extern "C" SEXP PairTracker__link(SEXP xptr, SEXP i, SEXP j, SEXP link)
{
    BEGIN_RCPP
    Rcpp::XPtr<PairTracker> ptr(xptr);
    return ptr->link(Rcpp::as<int>(i), Rcpp::as<int>(j), Rcpp::as<int>(link));
    END_RCPP
}

extern "C" SEXP PairTracker__getUnseenPair(SEXP xptr)
{
    BEGIN_RCPP
    Rcpp::XPtr<PairTracker> ptr(xptr);
    return ptr->getUnseenPair();
    END_RCPP
}

} // namespace dtwclust
