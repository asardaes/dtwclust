#ifndef DTWCLUST_H_
#define DTWCLUST_H_

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// see https://stackoverflow.com/questions/8551418/c-preprocessor-macro-for-returning-a-string-repeated-a-certain-number-of-times
#define DTWCLUST_REP0(X)
#define DTWCLUST_REP1(X) X
#define DTWCLUST_REP2(X) DTWCLUST_REP1(X), X
#define DTWCLUST_REP3(X) DTWCLUST_REP2(X), X
#define DTWCLUST_REP4(X) DTWCLUST_REP3(X), X
#define DTWCLUST_REP5(X) DTWCLUST_REP4(X), X
#define DTWCLUST_REP6(X) DTWCLUST_REP5(X), X
#define DTWCLUST_REP7(X) DTWCLUST_REP6(X), X
#define DTWCLUST_REP8(X) DTWCLUST_REP7(X), X
#define DTWCLUST_REP9(X) DTWCLUST_REP8(X), X
#define DTWCLUST_REP10(X) DTWCLUST_REP9(X), X

#define DTWCLUST_COMMA0
#define DTWCLUST_COMMA1 ,
#define DTWCLUST_COMMA2 ,
#define DTWCLUST_COMMA3 ,
#define DTWCLUST_COMMA4 ,
#define DTWCLUST_COMMA5 ,
#define DTWCLUST_COMMA6 ,
#define DTWCLUST_COMMA7 ,
#define DTWCLUST_COMMA8 ,
#define DTWCLUST_COMMA9 ,

// below will only work until 19 arguments
#define DTWCLUST_REP(TENS, ONES, X)\
    DTWCLUST_REP##TENS(DTWCLUST_REP10(X) DTWCLUST_COMMA##ONES) \
    DTWCLUST_REP##ONES(X)

#define DTWCLUST_GET_FUN(TENS, ONES, NAME) \
    static SEXP(*fun)(DTWCLUST_REP(TENS, ONES, SEXP)) = NULL; \
    if (fun == NULL) \
        fun = (SEXP(*)(DTWCLUST_REP(TENS, ONES, SEXP))) R_GetCCallable("dtwclust", NAME)

#ifdef __cplusplus
namespace dtwclust {

extern "C" {
#endif

SEXP dba(SEXP X, SEXP CENT,
         SEXP MAX_ITER, SEXP DELTA, SEXP TRACE,
         SEXP MV, SEXP MV_VER, SEXP DOTS, SEXP NUM_THREADS)
{
    DTWCLUST_GET_FUN(0,9, "dba");
    return fun(X, CENT, MAX_ITER, DELTA, TRACE, MV, MV_VER, DOTS, NUM_THREADS);
}

SEXP dtw_basic(SEXP x, SEXP y, SEXP window,
               SEXP m, SEXP n, SEXP num_var,
               SEXP norm, SEXP step, SEXP backtrack, SEXP normalize,
               SEXP distmat)
{
    DTWCLUST_GET_FUN(1,1, "dtw_basic");
    return fun(x, y, window, m, n, num_var, norm, step, backtrack, normalize, distmat);
}

SEXP dtw_lb(SEXP X, SEXP Y, SEXP D, SEXP MARGIN, SEXP DOTS, SEXP NUM_THREADS)
{
    DTWCLUST_GET_FUN(0,6, "dtw_lb");
    return fun(X, Y, D, MARGIN, DOTS, NUM_THREADS);
}

SEXP envelope(SEXP series, SEXP window)
{
    DTWCLUST_GET_FUN(0,2, "envelope");
    return fun(series, window);
}

SEXP lbi(SEXP X, SEXP Y, SEXP WINDOW, SEXP P, SEXP L, SEXP U)
{
    DTWCLUST_GET_FUN(0,6, "lbi");
    return fun(X, Y, WINDOW, P, L, U);
}

SEXP lbk(SEXP X, SEXP P, SEXP L, SEXP U)
{
    DTWCLUST_GET_FUN(0,4, "lbk");
    return fun(X, P, L, U);
}

SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP num_var, SEXP sigma, SEXP window, SEXP logs)
{
    DTWCLUST_GET_FUN(0,8, "logGAK");
    return fun(x, y, nx, ny, num_var, sigma, window, logs);
}

SEXP sdtw_cent(SEXP SERIES, SEXP CENTROID,
               SEXP GAMMA, SEXP WEIGHTS,
               SEXP MV, SEXP NUM_THREADS)
{
    DTWCLUST_GET_FUN(0,6, "sdtw_cent");
    return fun(SERIES, CENTROID, GAMMA, WEIGHTS, MV, NUM_THREADS);
}

SEXP soft_dtw(SEXP X, SEXP Y, SEXP GAMMA, SEXP COSTMAT, SEXP MV)
{
    DTWCLUST_GET_FUN(0,5, "soft_dtw");
    return fun(X, Y, GAMMA, COSTMAT, MV);
}

SEXP tadpole(SEXP X, SEXP K, SEXP DC, SEXP DTW_ARGS,
             SEXP LB, SEXP UB, SEXP TRACE,
             SEXP LIST, SEXP NUM_THREADS)
{
    DTWCLUST_GET_FUN(0,9, "tadpole");
    return fun(X, K, DC, DTW_ARGS, LB, UB, TRACE, LIST, NUM_THREADS);
}

#ifdef __cplusplus
} // extern C

} // namespace dtwclust
#endif

#endif // DTWCLUST_H_
