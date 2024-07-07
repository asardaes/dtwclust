roxygen_error_check_param <- function() {
    "Logical indicating whether the function should try to detect inconsistencies
    and give more informative errors messages. Also used internally to avoid repeating checks."
}

roxygen_window_details <- function() {
    "The windowing constraint uses a centered window.
The calculations expect a value in `window.size` that represents the distance between the point considered and one of the edges of the window.
Therefore, if, for example, `window.size = 10`, the warping for an observation \\eqn{x_i} considers the points between \\eqn{x_{i-10}} and \\eqn{x_{i+10}},
resulting in `10(2) + 1 = 21` observations falling within the window."
}

roxygen_proxy_section <- function() {
    "Proxy version:

The version registered with [proxy::dist()] is custom (`loop = FALSE` in [proxy::pr_DB]).
The custom function handles multi-threaded parallelization directly with [RcppParallel][RcppParallel::RcppParallel-package].
It uses all available threads by default (see [RcppParallel::defaultNumThreads()]),
but this can be changed by the user with [RcppParallel::setThreadOptions()].

An exception to the above is when it is called within a [`foreach`][foreach::foreach] parallel loop **made by dtwclust**.
If the parallel workers do not have the number of threads explicitly specified,
this function will default to 1 thread per worker.
See the parallelization vignette for more information - `browseVignettes(\"dtwclust\")`"
}

roxygen_proxy_symmetric <- function() {
    "It also includes symmetric optimizations to calculate only half a distance matrix when appropriate---only one list of series should be provided in `x`.
Starting with version 6.0.0, this optimization means that the function returns an array with the lower triangular values of the distance matrix,
similar to what [stats::dist()] does;
see [DistmatLowerTriangular-class] for a helper to access elements as it if were a normal matrix.
If you want to avoid this optimization, call [proxy::dist] by giving the same list of series in both `x` and `y`."
}

roxygen_rcpp_parallel_section <- function() {
    "Parallel Computing:
Please note that running tasks in parallel does **not** guarantee faster computations.
The overhead introduced is sometimes too large, and it's better to run tasks sequentially.

This function uses the [RcppParallel][RcppParallel::RcppParallel-package] package for parallelization.
It uses all available threads by default (see [RcppParallel::defaultNumThreads()]),
but this can be changed by the user with [RcppParallel::setThreadOptions()].

An exception to the above is when it is called within a [`foreach`][foreach::foreach] parallel loop **made by dtwclust**.
If the parallel workers do not have the number of threads explicitly specified,
this function will default to 1 thread per worker.
See the parallelization vignette for more information - `browseVignettes(\"dtwclust\")`"
}
