#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

void sbd_loop_pairwise(MatrixAccessor<double>& dist, const int fftlen,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                       const Rcpp::List& endpoints)
{
    arma::vec cc_seq_truncated(fftlen);
    int i_start = endpoints["start"], i_end = endpoints["end"];

    for (int i = i_start - 1; i < i_end; i++) {
        R_CheckUserInterrupt();

        // in two steps to avoid disambiguation
        Rcpp::NumericVector x_rcpp(X[i]);
        Rcpp::NumericVector y_rcpp(Y[i]);
        Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
        Rcpp::ComplexVector ffty_rcpp(FFTY[i]);
        arma::vec x(x_rcpp), y(y_rcpp);
        arma::cx_vec fftx(fftx_rcpp), ffty(ffty_rcpp);

        // already normalizes by length
        arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

        // reorder truncated sequence
        int id = 0;
        for (unsigned int j = fftlen - y.size() + 1; j < cc_seq.size(); j++) {
            cc_seq_truncated[id] = cc_seq[j];
            id++;
        }
        for (unsigned int j = 0; j < x.size(); j++) {
            cc_seq_truncated[id] = cc_seq[j];
            id++;
        }

        // get max
        double cc_max = R_NegInf;
        double den = arma::norm(x) * arma::norm(y);
        for (int j = 0; j < id; j++) {
            double this_cc = cc_seq_truncated[j] / den;
            if (this_cc > cc_max) cc_max = this_cc;
        }

        // bigmemory operator[][] is backwards
        dist[0][i] = 1 - cc_max;
    }
}

// =================================================================================================
/* symmetric case */
// =================================================================================================

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

            // already normalizes by length
            arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

            // reorder truncated sequence
            int id = 0;
            for (unsigned int k = fftlen - y.size() + 1; k < cc_seq.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }
            for (unsigned int k = 0; k < x.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }

            // get max
            double cc_max = R_NegInf;
            double den = arma::norm(x) * y_norm;
            for (int k = 0; k < id; k++) {
                double this_cc = cc_seq_truncated[k] / den;
                if (this_cc > cc_max) cc_max = this_cc;
            }

            // bigmemory operator[][] is backwards
            cc_max = 1 - cc_max;
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

void sbd_loop_general(MatrixAccessor<double>& dist, const int fftlen,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                      const Rcpp::List& endpoints)
{
    arma::vec cc_seq_truncated(fftlen);
    int j_start = endpoints["start"], j_end = endpoints["end"];

    for (int i = 0; i < X.length(); i++) {
        // in two steps to avoid disambiguation
        Rcpp::NumericVector x_rcpp(X[i]);
        Rcpp::ComplexVector fftx_rcpp(FFTX[i]);
        arma::vec x(x_rcpp);
        arma::cx_vec fftx(fftx_rcpp);
        double x_norm = arma::norm(x);

        for (int j = j_start - 1; j < j_end; j++) {
            R_CheckUserInterrupt();

            // in two steps to avoid disambiguation
            Rcpp::NumericVector y_rcpp(Y[j]);
            Rcpp::ComplexVector ffty_rcpp(FFTY[j]);
            arma::vec y(y_rcpp);
            arma::cx_vec ffty(ffty_rcpp);

            // already normalizes by length
            arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

            // reorder truncated sequence
            int id = 0;
            for (unsigned int k = fftlen - y.size() + 1; k < cc_seq.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }
            for (unsigned int k = 0; k < x.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }

            // get max
            double cc_max = R_NegInf;
            double den = x_norm * arma::norm(y);
            for (int k = 0; k < id; k++) {
                double this_cc = cc_seq_truncated[k] / den;
                if (this_cc > cc_max) cc_max = this_cc;
            }

            // bigmemory operator[][] is backwards
            dist[j][i] = 1 - cc_max;
        }
    }
}

// =================================================================================================
/* gateway function */
// =================================================================================================

RcppExport SEXP sbd_loop(SEXP D, SEXP X, SEXP Y, SEXP FFTX, SEXP FFTY,
                         SEXP FFTLEN, SEXP SYMMETRIC, SEXP PAIRWISE, SEXP ENDPOINTS)
{
    BEGIN_RCPP
    Rcpp::XPtr<BigMatrix> dist_ptr(D);
    MatrixAccessor<double> dist(*dist_ptr);

    if (Rcpp::as<bool>(PAIRWISE))
        sbd_loop_pairwise(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY, ENDPOINTS);
    else if (Rcpp::as<bool>(SYMMETRIC))
        sbd_loop_symmetric(dist, Rcpp::as<int>(FFTLEN), X, FFTX, FFTY, ENDPOINTS);
    else
        sbd_loop_general(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY, ENDPOINTS);

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
