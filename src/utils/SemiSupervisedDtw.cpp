#include "utils.h"

#include <cmath> // std::round

#include <RcppArmadillo.h>

#include "UndirectedGraph.h"

namespace dtwclust {

#define DONT_KNOW -1
#define CANNOT_LINK 0
#define MUST_LINK 1

class SemiSupervisedDtw
{
public:
    SemiSupervisedDtw(int max_size)
        : must_link_(max_size)
        , cannot_link_(max_size)
        , dont_know_(max_size)
        , max_size_(max_size)
    { }

    SEXP link(const int i, const int j, const int link_type) {
        if (i < 1 || i > max_size_ || j < 1 || j > max_size_) Rcpp::stop("Invalid indices provided");
        switch(link_type)
        {
        case DONT_KNOW:
            dont_know_.linkVertices(i,j);
            return Rcpp::wrap(dont_know_.isConnected());
        case CANNOT_LINK:
            cannot_link_.linkVertices(i,j);
            return Rcpp::wrap(cannot_link_.isConnected());
        case MUST_LINK:
            must_link_.linkVertices(i,j);
            return Rcpp::wrap(must_link_.isConnected());
        }
    }

    SEXP getUnseenPair() {
        if (!(this->optionsAvailable())) return R_NilValue;
        Rcpp::IntegerVector pair(2);
        Rcpp::RNGScope rng_scope;
        bool seen = true;
        while (seen) {
            Rcpp::checkUserInterrupt();
            seen = false;
            pair[0] = std::round(R::runif(1, max_size_));
            pair[1] = std::round(R::runif(1, max_size_));
            while (pair[0] == pair[1]) pair[1] = std::round(R::runif(1, max_size_));
            if (must_link_.areNeighbors(pair[0], pair[1])) seen = true;
            if (cannot_link_.areNeighbors(pair[0], pair[1])) seen = true;
            if (dont_know_.areNeighbors(pair[0], pair[1])) seen = true;
        }
        return pair;
    }

private:
    UndirectedGraph must_link_, cannot_link_, dont_know_;
    int max_size_;

    bool optionsAvailable() {
        if (must_link_.isConnected() || cannot_link_.isConnected() || dont_know_.isConnected())
            return false;

        int total_edges = must_link_.numEdges() + cannot_link_.numEdges() + dont_know_.numEdges();
        if (total_edges >= max_size_ * (max_size_ - 1) / 2)
            return false;

        return true;
    }
};

RcppExport SEXP SemiSupervisedDtw__new(SEXP max_size)
{
    BEGIN_RCPP
    Rcpp::XPtr<SemiSupervisedDtw> ptr(new SemiSupervisedDtw(Rcpp::as<int>(max_size)), true);
    return ptr;
    END_RCPP
}

RcppExport SEXP SemiSupervisedDtw__link(SEXP xptr, SEXP i, SEXP j, SEXP link)
{
    BEGIN_RCPP
    Rcpp::XPtr<SemiSupervisedDtw> ptr(xptr);
    return ptr->link(Rcpp::as<int>(i), Rcpp::as<int>(j), Rcpp::as<int>(link));
    END_RCPP
}

RcppExport SEXP SemiSupervisedDtw__getUnseenPair(SEXP xptr)
{
    BEGIN_RCPP
    Rcpp::XPtr<SemiSupervisedDtw> ptr(xptr);
    return ptr->getUnseenPair();
    END_RCPP
}

} // namespace dtwclust
