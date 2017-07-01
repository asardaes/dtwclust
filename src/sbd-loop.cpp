#include <RcppArmadillo.h>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

void sbd_loop_pairwise(Rcpp::NumericVector& dist, const int len,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& FFTX, const Rcpp::List& FFTY)
{
    arma::vec cc_seq_truncated(len);
    for (int i = 0; i < X.length(); i++) {
        arma::vec x(Rcpp::NumericVector(X[i]));
        arma::vec y(Rcpp::NumericVector(Y[i]));
        arma::cx_vec fftx(Rcpp::ComplexVector(FFTX[i]));
        arma::cx_vec ffty(Rcpp::ComplexVector(FFTY[i]));

        // already normalizes by length
        arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

        // reorder truncated sequence
        int id = 0;
        for (int j = len - y.size() + 1; j < cc_seq.size(); j++) {
            cc_seq_truncated[id] = cc_seq[j];
            id++;
        }
        for (int j = 0; j < x.size(); j++) {
            cc_seq_truncated[id] = cc_seq[j];
            id++;
        }
        while (id < len) {
            cc_seq_truncated[id] = 0;
            id++;
        }

        // get max
        double cc_max = R_NegInf;
        double den = arma::norm(x) * arma::norm(y);
        for (int j = 0; j < len; j++) {
            double this_cc = cc_seq_truncated[j] / den;
            if (this_cc > cc_max) cc_max = this_cc;
        }

        dist[i] = 1 - cc_max;
    }
}

// =================================================================================================
/* general case */
// =================================================================================================

void sbd_loop_general(Rcpp::NumericMatrix& dist, const int len,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& FFTX, const Rcpp::List& FFTY)
{
    arma::vec cc_seq_truncated(len);
    for (int i = 0; i < X.length(); i++) {
        for (int j = 0; j < Y.length(); j++) {
            arma::vec x(Rcpp::NumericVector(X[i]));
            arma::vec y(Rcpp::NumericVector(Y[j]));
            arma::cx_vec fftx(Rcpp::ComplexVector(FFTX[i]));
            arma::cx_vec ffty(Rcpp::ComplexVector(FFTY[j]));

            // already normalizes by length
            arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

            // reorder truncated sequence
            int id = 0;
            for (int k = len - y.size() + 1; k < cc_seq.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }
            for (int k = 0; k < x.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }
            while (id < len) {
                cc_seq_truncated[id] = 0;
                id++;
            }

            // get max
            double cc_max = R_NegInf;
            double den = arma::norm(x) * arma::norm(y);
            for (int k = 0; k < len; k++) {
                double this_cc = cc_seq_truncated[k] / den;
                if (this_cc > cc_max) cc_max = this_cc;
            }

            dist(i,j) = 1 - cc_max;
        }
    }
}

// =================================================================================================
/* main gateway function */
// =================================================================================================

RcppExport SEXP sbd_loop(SEXP D, SEXP X, SEXP Y, SEXP FFTX, SEXP FFTY, SEXP FFTLEN, SEXP PAIRWISE)
{
    BEGIN_RCPP
    if (Rcpp::as<bool>(PAIRWISE)) {
        Rcpp::NumericVector dist(D);
        sbd_loop_pairwise(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY);

    } else {
        Rcpp::NumericMatrix dist(D);
        sbd_loop_general(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY);
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
