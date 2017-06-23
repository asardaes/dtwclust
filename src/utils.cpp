#include <Rcpp.h>
#include "dtwclustpp.h"
#include "dtwclust.h"

namespace dtwclust {

// =================================================================================================
/* dtw_basic proxy */
// =================================================================================================

double dtwb(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& dots)
{
    SEXP X = PROTECT(Rcpp::wrap(x));
    SEXP Y = PROTECT(Rcpp::wrap(y));
    SEXP NX = PROTECT(Rcpp::wrap(x.length()));
    SEXP NY = PROTECT(Rcpp::wrap(y.length()));
    SEXP NV = PROTECT(Rcpp::wrap(1));

    SEXP window = dots["window.size"];
    SEXP norm = dots["norm"];
    SEXP step = dots["step.pattern"];
    SEXP backtrack = dots["backtrack"];
    SEXP gcm = dots["gcm"];

    double d = Rcpp::as<double>(dtw_basic(X, Y, window, NX, NY, NV, norm, step, backtrack, gcm));
    UNPROTECT(5);
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
