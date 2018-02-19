#include "utils.h"

#include <cmath> // std::round

#include <RcppArmadillo.h>

#include "UndirectedGraph.h"

namespace dtwclust {

class SemiSupervisedDtw
{
public:
    SemiSupervisedDtw(int max_size)
        : graph_(max_size)
        , max_size_(max_size)
    { }

    SEXP link(const int i, const int j, const int link_type) {
        if (i < 1 || j < 1) Rcpp::stop("Invalid indices provided");
        graph_.linkVertices(i, j, link_type);
        return Rcpp::wrap(graph_.isConnected());
    }

    SEXP getUnseenPair() {
        if (graph_.isConnected()) return R_NilValue;
        Rcpp::IntegerVector pair(2);
        Rcpp::RNGScope rng_scope;
        bool seen = true;
        while (seen) {
            Rcpp::checkUserInterrupt();
            seen = false;
            pair[0] = std::round(R::runif(1, max_size_));
            pair[1] = std::round(R::runif(1, max_size_));
            while (pair[0] == pair[1]) pair[1] = std::round(R::runif(1, max_size_));
            if (graph_.areNeighbors(pair[0], pair[1])) seen = true;
        }
        return pair;
    }

private:
    UndirectedGraph graph_;
    int max_size_;
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
