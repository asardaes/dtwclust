#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclust++.h"

namespace dtwclust {

// =================================================================================================
/* calculate the shape-based distance */
// =================================================================================================

double sbd_core(const arma::cx_vec& fftx, const arma::cx_vec& ffty,
                const arma::vec& x, const arma::vec& y,
                const int fftlen, const double x_norm, const double y_norm,
                arma::vec& cc_seq_truncated)
{
    // already normalizes by length
    arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

    // reorder truncated sequence
    int id = 0;
    for (unsigned int i = fftlen - y.size() + 1; i < cc_seq.size(); i++) {
        cc_seq_truncated[id] = cc_seq[i];
        id++;
    }
    for (unsigned int i = 0; i < x.size(); i++) {
        cc_seq_truncated[id] = cc_seq[i];
        id++;
    }

    // get max
    double cc_max = R_NegInf;
    double den = x_norm * y_norm;
    for (int i = 0; i < id; i++) {
        double this_cc = cc_seq_truncated[i] / den;
        if (this_cc > cc_max) cc_max = this_cc;
    }

    return 1 - cc_max;
}

// =================================================================================================
/* pairwise case */
// =================================================================================================

// R matrix
void sbd_loop_pairwise(Rcpp::NumericMatrix& dist, const int fftlen,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                       int index)
{
    arma::vec cc_seq_truncated(fftlen);
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();

        // in two steps to avoid disambiguation
        Rcpp::NumericVector x_rcpp(X[i]);
        Rcpp::NumericVector y_rcpp(Y[i]);
        Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[i]);
        arma::vec x(x_rcpp), y(y_rcpp);
        arma::cx_vec fftx(fftx_rcpp), ffty(ffty_rcpp);

        double x_norm = arma::norm(x);
        double y_norm = arma::norm(y);
        dist(index++, 0) = sbd_core(fftx, ffty, x, y, fftlen, x_norm, y_norm, cc_seq_truncated);
    }
}

// big.matrix
void sbd_loop_pairwise(MatrixAccessor<double>& dist, const int fftlen,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                       int index)
{
    arma::vec cc_seq_truncated(fftlen);
    index--;
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();

        // in two steps to avoid disambiguation
        Rcpp::NumericVector x_rcpp(X[i]);
        Rcpp::NumericVector y_rcpp(Y[i]);
        Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[i]);
        arma::vec x(x_rcpp), y(y_rcpp);
        arma::cx_vec fftx(fftx_rcpp), ffty(ffty_rcpp);

        // bigmemory operator[][] is backwards
        double x_norm = arma::norm(x);
        double y_norm = arma::norm(y);
        dist[0][index++] = sbd_core(fftx, ffty, x, y, fftlen, x_norm, y_norm, cc_seq_truncated);
    }
}

// =================================================================================================
/* symmetric case */
// =================================================================================================

// R matrix
void sbd_loop_symmetric(Rcpp::NumericMatrix& dist, const int fftlen,
                        const Rcpp::List& X, const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                        const Rcpp::List& endpoints)
{
    arma::vec cc_seq_truncated(fftlen);
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        // in two steps to avoid disambiguation
        Rcpp::NumericVector y_rcpp(X[j]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[j]);
        arma::vec y(y_rcpp);
        arma::cx_vec ffty(ffty_rcpp);
        double y_norm = arma::norm(y);

        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = X.length();

        while (i < i_max) {
            R_CheckUserInterrupt();

            // in two steps to avoid disambiguation
            Rcpp::NumericVector x_rcpp(X[i]);
            Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
            arma::vec x(x_rcpp);
            arma::cx_vec fftx(fftx_rcpp);

            double x_norm = arma::norm(x);
            double cc_max = sbd_core(fftx, ffty, x, y, fftlen, x_norm, y_norm, cc_seq_truncated);

            dist(i,j) = cc_max;
            dist(j,i) = cc_max;
            i++;
        }
        j++;
        i = j + 1;
    }
}

// big.matrix
void sbd_loop_symmetric(MatrixAccessor<double>& dist, const int fftlen,
                        const Rcpp::List& X, const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                        const Rcpp::List& endpoints)
{
    arma::vec cc_seq_truncated(fftlen);
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        // in two steps to avoid disambiguation
        Rcpp::NumericVector y_rcpp(X[j]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[j]);
        arma::vec y(y_rcpp);
        arma::cx_vec ffty(ffty_rcpp);
        double y_norm = arma::norm(y);

        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = X.length();

        while (i < i_max) {
            R_CheckUserInterrupt();

            // in two steps to avoid disambiguation
            Rcpp::NumericVector x_rcpp(X[i]);
            Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
            arma::vec x(x_rcpp);
            arma::cx_vec fftx(fftx_rcpp);

            double x_norm = arma::norm(x);
            double cc_max = sbd_core(fftx, ffty, x, y, fftlen, x_norm, y_norm, cc_seq_truncated);

            // bigmemory operator[][] is backwards
            dist[i][j] = cc_max;
            dist[j][i] = cc_max;
            i++;
        }
        j++;
        i = j + 1;
    }
}

// =================================================================================================
/* general case */
// =================================================================================================

// R matrix
void sbd_loop_general(Rcpp::NumericMatrix& dist, const int fftlen,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                      int index)
{
    arma::vec cc_seq_truncated(fftlen);
    index--;
    for (int j = 0; j < Y.length(); j++) {
        // in two steps to avoid disambiguation
        Rcpp::NumericVector y_rcpp(Y[j]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[j]);
        arma::vec y(y_rcpp);
        arma::cx_vec ffty(ffty_rcpp);
        double y_norm = arma::norm(y);

        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();

            // in two steps to avoid disambiguation
            Rcpp::NumericVector x_rcpp(X[i]);
            Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
            arma::vec x(x_rcpp);
            arma::cx_vec fftx(fftx_rcpp);

            double x_norm = arma::norm(x);

            dist(i,index) = sbd_core(fftx, ffty, x, y, fftlen, x_norm, y_norm, cc_seq_truncated);
        }
        index++;
    }
}

// big.matrix
void sbd_loop_general(MatrixAccessor<double>& dist, const int fftlen,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                      int index)
{
    arma::vec cc_seq_truncated(fftlen);
    index--;
    for (int j = 0; j < Y.length(); j++) {
        // in two steps to avoid disambiguation
        Rcpp::NumericVector y_rcpp(Y[j]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[j]);
        arma::vec y(y_rcpp);
        arma::cx_vec ffty(ffty_rcpp);
        double y_norm = arma::norm(y);

        for (int i = 0; i < X.length(); i++) {
            R_CheckUserInterrupt();

            // in two steps to avoid disambiguation
            Rcpp::NumericVector x_rcpp(X[i]);
            Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
            arma::vec x(x_rcpp);
            arma::cx_vec fftx(fftx_rcpp);

            double x_norm = arma::norm(x);

            // bigmemory operator[][] is backwards
            dist[index][i] = sbd_core(fftx, ffty, x, y, fftlen, x_norm, y_norm, cc_seq_truncated);
        }
        index++;
    }
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP sbd_loop(SEXP D, SEXP X, SEXP Y, SEXP FFTX, SEXP FFTY,
                         SEXP FFTLEN, SEXP SYMMETRIC, SEXP PAIRWISE, SEXP ENDPOINTS, SEXP BIGMAT)
{
    BEGIN_RCPP
    if (Rcpp::as<bool>(BIGMAT)) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        MatrixAccessor<double> dist(*dist_ptr);

        if (Rcpp::as<bool>(PAIRWISE))
            sbd_loop_pairwise(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY, Rcpp::as<int>(ENDPOINTS));
        else if (Rcpp::as<bool>(SYMMETRIC))
            sbd_loop_symmetric(dist, Rcpp::as<int>(FFTLEN), X, FFTX, FFTY, ENDPOINTS);
        else
            sbd_loop_general(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY, Rcpp::as<int>(ENDPOINTS));

    } else {
        Rcpp::NumericMatrix dist(D);

        if (Rcpp::as<bool>(PAIRWISE))
            sbd_loop_pairwise(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY, Rcpp::as<int>(ENDPOINTS));
        else if (Rcpp::as<bool>(SYMMETRIC))
            sbd_loop_symmetric(dist, Rcpp::as<int>(FFTLEN), X, FFTX, FFTY, ENDPOINTS);
        else
            sbd_loop_general(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY, Rcpp::as<int>(ENDPOINTS));
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
