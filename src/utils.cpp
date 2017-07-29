#include <Rcpp.h>
#include "dtwclust++.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* gak proxy */
// =================================================================================================

double gak(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& dots)
{
    SEXP NX = PROTECT(Rcpp::wrap(x.length()));
    SEXP NY = PROTECT(Rcpp::wrap(y.length()));
    SEXP NV = PROTECT(Rcpp::wrap(1));

    SEXP sigma = dots["sigma"];
    SEXP window = dots["window.size"];
    SEXP logs = dots["logs"];

    double d = Rcpp::as<double>(logGAK(x, y, NX, NY, NV, sigma, window, logs));
    UNPROTECT(3);
    return d;
}

double gak(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y, const Rcpp::List& dots)
{
    SEXP NX = PROTECT(Rcpp::wrap(x.nrow()));
    SEXP NY = PROTECT(Rcpp::wrap(y.nrow()));
    SEXP NV = PROTECT(Rcpp::wrap(x.ncol()));

    SEXP sigma = dots["sigma"];
    SEXP window = dots["window.size"];
    SEXP logs = dots["logs"];

    double d = Rcpp::as<double>(logGAK(x, y, NX, NY, NV, sigma, window, logs));
    UNPROTECT(3);
    return d;
}

// =================================================================================================
/* dtw_basic proxy */
// =================================================================================================

double dtwb(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& dots)
{
    SEXP NX = PROTECT(Rcpp::wrap(x.length()));
    SEXP NY = PROTECT(Rcpp::wrap(y.length()));
    SEXP NV = PROTECT(Rcpp::wrap(1));

    SEXP window = dots["window.size"];
    SEXP norm = dots["norm"];
    SEXP step = dots["step.pattern"];
    SEXP backtrack = dots["backtrack"];
    SEXP gcm = dots["gcm"];

    double d = Rcpp::as<double>(dtw_basic(x, y, window, NX, NY, NV, norm, step, backtrack, gcm));
    UNPROTECT(3);
    return d;
}

double dtwb(const Rcpp::NumericMatrix& x, const Rcpp::NumericMatrix& y, const Rcpp::List& dots)
{
    SEXP NX = PROTECT(Rcpp::wrap(x.nrow()));
    SEXP NY = PROTECT(Rcpp::wrap(y.nrow()));
    SEXP NV = PROTECT(Rcpp::wrap(x.ncol()));

    SEXP window = dots["window.size"];
    SEXP norm = dots["norm"];
    SEXP step = dots["step.pattern"];
    SEXP backtrack = dots["backtrack"];
    SEXP gcm = dots["gcm"];

    double d = Rcpp::as<double>(dtw_basic(x, y, window, NX, NY, NV, norm, step, backtrack, gcm));
    UNPROTECT(3);
    return d;
}

// =================================================================================================
/* for Rcpp::Rcout */
// =================================================================================================

void Rflush()
{
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
}

} // namespace dtwclust
