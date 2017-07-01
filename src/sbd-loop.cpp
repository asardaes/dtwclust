#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "dtwclustpp.h"

namespace dtwclust {

// =================================================================================================
/* pairwise case */
// =================================================================================================

void sbd_loop_pairwise(Rcpp::NumericVector& dist, const int fftlen,
                       const Rcpp::List& X, const Rcpp::List& Y,
                       const Rcpp::List& FFTX, const Rcpp::List& FFTY)
{
    arma::vec cc_seq_truncated(fftlen);
    for (int i = 0; i < X.length(); i++) {
        R_CheckUserInterrupt();
        arma::vec x(Rcpp::NumericVector(X[i]));
        arma::vec y(Rcpp::NumericVector(Y[i]));
        arma::cx_vec fftx(Rcpp::ComplexVector(FFTX[i]));
        arma::cx_vec ffty(Rcpp::ComplexVector(FFTY[i]));

        // already normalizes by length
        arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

        // reorder truncated sequence
        int id = 0;
        for (int j = fftlen - y.size() + 1; j < cc_seq.size(); j++) {
            cc_seq_truncated[id] = cc_seq[j];
            id++;
        }
        for (int j = 0; j < x.size(); j++) {
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

        dist[i] = 1 - cc_max;
    }
}

// =================================================================================================
/* symmetric case */
// =================================================================================================

void sbd_loop_symmetric(const Rcpp::XPtr<BigMatrix>& dist_ptr, const int fftlen,
                        const Rcpp::List& X, const Rcpp::List& FFTX, const Rcpp::List& FFTY,
                        const Rcpp::List& endpoints)
{
    MatrixAccessor<double> dist(*dist_ptr);
    arma::vec cc_seq_truncated(fftlen);
    Rcpp::List start = Rcpp::as<Rcpp::List>(endpoints["start"]);
    Rcpp::List end = Rcpp::as<Rcpp::List>(endpoints["end"]);
    int i_start = start["i"], i_end = end["i"];
    int j_start = start["j"], j_end = end["j"];

    int i = i_start - 1, j = j_start - 1;
    while (j < j_end) {
        int i_max;
        if (j == (j_end - 1))
            i_max = i_end;
        else
            i_max = X.length();

        while (i < i_max) {
            R_CheckUserInterrupt();
            arma::vec x(Rcpp::NumericVector(X[i]));
            arma::vec y(Rcpp::NumericVector(X[j]));
            arma::cx_vec fftx(Rcpp::ComplexVector(FFTX[i]));
            arma::cx_vec ffty(Rcpp::ComplexVector(FFTY[j]));

            // already normalizes by length
            arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

            // reorder truncated sequence
            int id = 0;
            for (int k = fftlen - y.size() + 1; k < cc_seq.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }
            for (int k = 0; k < x.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }

            // get max
            double cc_max = R_NegInf;
            double den = arma::norm(x) * arma::norm(y);
            for (int k = 0; k < id; k++) {
                double this_cc = cc_seq_truncated[k] / den;
                if (this_cc > cc_max) cc_max = this_cc;
            }

            // bigmemory operator[][] is backwards, but it doesn't matter
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

void sbd_loop_general(Rcpp::NumericMatrix& dist, const int fftlen,
                      const Rcpp::List& X, const Rcpp::List& Y,
                      const Rcpp::List& FFTX, const Rcpp::List& FFTY)
{
    arma::vec cc_seq_truncated(fftlen);
    for (int i = 0; i < X.length(); i++) {
        for (int j = 0; j < Y.length(); j++) {
            R_CheckUserInterrupt();
            arma::vec x(Rcpp::NumericVector(X[i]));
            arma::vec y(Rcpp::NumericVector(Y[j]));
            arma::cx_vec fftx(Rcpp::ComplexVector(FFTX[i]));
            arma::cx_vec ffty(Rcpp::ComplexVector(FFTY[j]));

            // already normalizes by length
            arma::vec cc_seq = arma::real(arma::ifft(fftx % ffty));

            // reorder truncated sequence
            int id = 0;
            for (int k = fftlen - y.size() + 1; k < cc_seq.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }
            for (int k = 0; k < x.size(); k++) {
                cc_seq_truncated[id] = cc_seq[k];
                id++;
            }

            // get max
            double cc_max = R_NegInf;
            double den = arma::norm(x) * arma::norm(y);
            for (int k = 0; k < id; k++) {
                double this_cc = cc_seq_truncated[k] / den;
                if (this_cc > cc_max) cc_max = this_cc;
            }

            dist(i,j) = 1 - cc_max;
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
    if (Rcpp::as<bool>(PAIRWISE)) {
        Rcpp::NumericVector dist(D);
        sbd_loop_pairwise(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY);

    } else if (Rcpp::as<bool>(SYMMETRIC)) {
        Rcpp::XPtr<BigMatrix> dist_ptr(D);
        sbd_loop_symmetric(dist_ptr, Rcpp::as<int>(FFTLEN), X, FFTX, FFTY, ENDPOINTS);

    } else {
        Rcpp::NumericMatrix dist(D);
        sbd_loop_general(dist, Rcpp::as<int>(FFTLEN), X, Y, FFTX, FFTY);
    }

    return R_NilValue;
    END_RCPP
}

} // namespace dtwclust
