# ==================================================================================================
# distance experiments using single series (1 to 1)
# ==================================================================================================

existing_objects <- ls(all.names = TRUE)

# --------------------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------------------

#' The dataset has series with different lengths. This threshold specifies the difference there
#' should be between series in the subset that is taken.
length_diff_threshold <- 20L

# Take only one variable of each character
series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
# Take only two samples of each character
series <- lapply(series, function(s) { s[1L:2L] })
# Get length of each character's time series
len <- sapply(series, function(s) { lengths(s)[1L] })
# Sort lengths
id_ascending <- sort(len, decreasing = FALSE, index.return = TRUE)

# Identify those characters that have length differences greater than the threshold specified
len <- id_ascending$x
id_ascending <- id_ascending$ix
dlen <- diff(len) < length_diff_threshold
while (any(dlen)) {
    rem <- which(dlen)[1L] + 1L
    len <- len[-rem]
    id_ascending <- id_ascending[-rem]
    dlen <- diff(len) < length_diff_threshold
}

# Get the resulting characters
series <- series[id_ascending]
# Also get the multivariate versions
series_mv <- lapply(multivariate_series, function(s) { s[1L:2L] })
series_mv <- series_mv[id_ascending]
# Window sizes for the experiments (where applicable)
window_sizes <- seq(from = 10L, to = 100L, by = 10L)
# Number of times each experiment will be repeated (by microbenchmark)
times <- if (short_experiments) 100L else 500L

#' NOTE: all single experiments are pretty much equivalent, only the first one is commented. They
#' are run within a new environment so that they don't change variables in the global environment
#' (this one).

# For metadata
t1 <- proc.time()

# --------------------------------------------------------------------------------------------------
# lb_keogh
# --------------------------------------------------------------------------------------------------

mycat("\tRunning lb_keogh experiments for single series\n")
dist_lbk_single <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series, function(this_series) {
            # Build expressions to evaluate, substituting window size
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    lb_keogh(x, y, .(window.size), error.check = FALSE)
                )
            })

            # Extract sample series
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            # Return data frame with results
            data.frame(distance = "lb_keogh",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median)
        })

        # Bind results for all series and return to global environment
        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# lb_improved
# --------------------------------------------------------------------------------------------------

mycat("\tRunning lb_improved experiments for single series\n")
dist_lbi_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    lb_improved(x, y, .(window.size), error.check = FALSE)
                )
            })

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "lb_improved",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# sbd
# --------------------------------------------------------------------------------------------------

mycat("\tRunning sbd experiments for single series\n")
dist_sbd_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(
                SBD(x, y, error.check = FALSE, return.shifted = FALSE),
                times = times, unit = "us"
            ))

            data.frame(distance = "sbd",
                       series_length = NROW(x),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# dtw univariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dtw experiments for single univariate series\n")
dist_dtw_univariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    dtw_basic(x, y, .(window.size), error.check = FALSE)
                )
            })

            expressions <- c(expressions, list(bquote(dtw_basic(x, y, error.check = FALSE))))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "dtw_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# dtw_multivariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dtw experiments for single multivariate series\n")
dist_dtw_multivariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series_mv, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    dtw_basic(x, y, .(window.size), error.check = FALSE)
                )
            })

            expressions <- c(expressions, list(bquote(dtw_basic(x, y, error.check = FALSE))))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "dtw_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# sdtw univariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning sdtw experiments for single univariate series\n")
dist_sdtw_univariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(
                sdtw(x, y, error.check = FALSE),
                times = times, unit = "us"
            ))

            data.frame(distance = "sdtw_univariate",
                       series_length = NROW(x),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# sdtw multivariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning sdtw experiments for single multivariate series\n")
dist_sdtw_multivariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series_mv, function(this_series) {
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(
                sdtw(x, y, error.check = FALSE),
                times = times, unit = "us"
            ))

            data.frame(distance = "sdtw_multivariate",
                       series_length = NROW(x),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# unnormalized gak univariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning unnormalized_gak experiments for single univariate series\n")
dist_unnormalized_gak_univariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    GAK(x, y, .(window.size), sigma = 100, normalize = FALSE, error.check = FALSE)
                )
            })

            expressions <- c(expressions, bquote(
                GAK(x, y, sigma = 100, normalize = FALSE, error.check = FALSE)
            ))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "unnormalized_gak_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# unnormalized gak multivariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning unnormalized_gak experiments for single multivariate series\n")
dist_unnormalized_gak_multivariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series_mv, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    GAK(x, y, .(window.size), sigma = 100, normalize = FALSE, error.check = FALSE)
                )
            })

            expressions <- c(expressions, bquote(
                GAK(x, y, sigma = 100, normalize = FALSE, error.check = FALSE)
            ))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "unnormalized_gak_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# normalized gak univariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning normalized_gak experiments for single univariate series\n")
dist_normalized_gak_univariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    GAK(x, y, .(window.size), sigma = 100, error.check = FALSE)
                )
            })

            expressions <- c(expressions, bquote(
                GAK(x, y, sigma = 100, error.check = FALSE)
            ))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "normalized_gak_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# normalized gak multivariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning normalized_gak experiments for single multivariate series\n")
dist_normalized_gak_multivariate_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series_mv, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    GAK(x, y, .(window.size), sigma = 100, error.check = FALSE)
                )
            })

            expressions <- c(expressions, bquote(
                GAK(x, y, sigma = 100, error.check = FALSE)
            ))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "normalized_gak_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median)
        })

        dplyr::bind_rows(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# aggregate
# --------------------------------------------------------------------------------------------------

dist_single_results <- dplyr::bind_rows(
    dist_lbk_single,
    dist_lbi_single,
    dist_sbd_single,
    dist_dtw_univariate_single,
    dist_dtw_multivariate_single,
    dist_sdtw_univariate_single,
    dist_sdtw_multivariate_single,
    dist_unnormalized_gak_univariate_single,
    dist_unnormalized_gak_multivariate_single,
    dist_normalized_gak_univariate_single,
    dist_normalized_gak_multivariate_single
)

# Make factor with the given order
dist_single_results$distance <- factor(dist_single_results$distance,
                                       levels = unique(dist_single_results$distance))

# Add some metadata
attr(dist_single_results, "proctime") <- proc.time() - t1
attr(dist_single_results, "times") <- times

# Clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "dist_single_results")))

# ==================================================================================================
# distance experiments using multiple series
# ==================================================================================================

existing_objects <- ls(all.names = TRUE)
t1 <- proc.time()

# --------------------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------------------

#' NOTE: these cases are almost the same as above, except we take all series for each character
#' (initially). Some new parameters are explained below.

length_diff_threshold <- 40L

series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
len <- sapply(series, function(s) { lengths(s)[1L] })
id_ascending <- sort(len, decreasing = FALSE, index.return = TRUE)

len <- id_ascending$x
id_ascending <- id_ascending$ix
dlen <- diff(len) < length_diff_threshold
while (any(dlen)) {
    rem <- which(dlen)[1L] + 1L
    len <- len[-rem]
    id_ascending <- id_ascending[-rem]
    dlen <- diff(len) < length_diff_threshold
}

series <- series[id_ascending]
series_mv <- multivariate_series[id_ascending]
# Window sizes for experiments that test more than one value
window_sizes <- seq(from = 30L, to = 90L, by = 30L)
# Window size for experiments that set a fixed value
window_size <- 50L

#' 'id_series' will have two columns, each one specifying the number of rows and columns the cross-
#' distance matrix should have. The short experiments only get square matrices.
if (short_experiments) {
    series <- series[1L:2L]
    series_mv <- series_mv[1L:2L]
    # Number of evaluations for each expression
    times <- 5L
    # Number of parallel threads to test
    num_threads_to_test <- c(4L)
    id_series <- cbind(seq(from = 10L, to = 100L, by = 10L),
                       seq(from = 10L, to = 100L, by = 10L))
} else {
    # Number of evaluations for each expression
    times <- 30L
    # Number of parallel threads to test
    num_threads_to_test <- c(1L, 2L, 4L)
    id_series <- rbind(
        expand.grid(seq(from = 20L, to = 100L, by = 20L), 20L),
        expand.grid(100L, seq(from = 40L, to = 100L, by = 20L)),
        cbind(Var1 = seq(from = 40L, to = 80L, by = 20L),
              Var2 = seq(from = 40L, to = 80L, by = 20L))
    )
    id_series <- id_series[order(id_series[,1L] * id_series[,2L]),]
}

cat("\n")

#' NOTE: the next experiments are also pretty much equivalent to each other, except some look along
#' different window sizes. Only the first one is commented. Also note that the proxy::dist version
#' of GAK is always normalized, that's why there aren't experiments for the unnormalized version.

# --------------------------------------------------------------------------------------------------
# lb_keogh
# --------------------------------------------------------------------------------------------------

mycat("\tRunning lb_keogh experiments for multiple series\n")
# Loop along number of parallel workers
dist_lbk_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")

    # Set number of threads to use
    RcppParallel::setThreadOptions(num_threads)

    # Loop along series
    benchmarks <- lapply(series, function(this_series) {
        # Build expressions for proxy::dist
        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "lbk",
                            window.size = .(window_size),
                            error.check = FALSE)
            )
        })

        # Evaluate expressions
        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        # Return data frame with results
        data.frame(distance = "lb_keogh",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = window_size,
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    # Bind results for all series and return
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# lb_improved
# --------------------------------------------------------------------------------------------------

mycat("\tRunning lb_improved experiments for multiple series\n")
dist_lbi_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    benchmarks <- lapply(series, function(this_series) {
        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "lbi",
                            window.size = .(window_size),
                            error.check = FALSE)
            )
        })

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "lb_improved",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = window_size,
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# dtw_lb
# --------------------------------------------------------------------------------------------------

#' NOTE: dtw_lb's experiments make more sense if the series in x and y are different, that's why
#' id_series is different here.

mycat("\tRunning dtw_lb experiments for multiple series\n")
dist_dtwlb_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    id_series <- rbind(
        expand.grid(seq(from = 10L, to = 50L, by = 10L), 10L),
        expand.grid(50L, seq(from = 20L, to = 100L, by = 10L)),
        cbind(Var1 = seq(from = 20L, to = 40L, by = 10L),
              Var2 = seq(from = 20L, to = 40L, by = 10L))
    )
    id_series <- id_series[order(id_series[,1L] * id_series[,2L]),]

    benchmarks <- lapply(series, function(this_series) {
        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[51L:.(50L + id_series[i, 2L])],
                            method = "dtw_lb",
                            window.size = .(window_size),
                            error.check = FALSE)
            )
        })

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "dtw_lb",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = window_size,
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# sbd
# --------------------------------------------------------------------------------------------------

mycat("\tRunning sbd experiments for multiple series\n")
dist_sbd_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    benchmarks <- lapply(series, function(this_series) {
        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "sbd",
                            error.check = FALSE)
            )
        })

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "sbd",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# dtw univariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dtw experiments for multiple univariate series\n")
dist_dtw_univariate_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    benchmarks <- lapply(series, function(this_series) {
        expressions <- lapply(window_sizes, function(window_size) {
            lapply(1L:nrow(id_series), function(i) {
                bquote(
                    proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                                y = this_series[1L:.(id_series[i, 2L])],
                                method = "dtw_basic",
                                window.size = .(window_size),
                                error.check = FALSE)
                )
            })
        })
        expressions <- unlist(expressions, recursive = FALSE)

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "dtw_univariate",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = rep(window_sizes, each = nrow(id_series)),
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# dtw multivariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning dtw experiments for multiple multivariate series\n")
dist_dtw_multivariate_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    benchmarks <- lapply(series_mv, function(this_series) {
        expressions <- lapply(window_sizes, function(window_size) {
            lapply(1L:nrow(id_series), function(i) {
                bquote(
                    proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                                y = this_series[1L:.(id_series[i, 2L])],
                                method = "dtw_basic",
                                window.size = .(window_size),
                                error.check = FALSE)
                )
            })
        })
        expressions <- unlist(expressions, recursive = FALSE)

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "dtw_multivariate",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = rep(window_sizes, each = nrow(id_series)),
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# sdtw univariate
# --------------------------------------------------------------------------------------------------

#' NOTE: GAK and soft-DTW are much more time consuming, so I only test univariate, less window
#' sizes, and less series. Based on the single experiments, this should be enough
#' to give an idea.

window_sizes <- c(20L, 40L)
id_series <- rbind(
    expand.grid(seq(from = 20L, to = 60L, by = 20L), 20L),
    expand.grid(60L, seq(from = 40L, to = 60L, by = 20L)),
    cbind(Var1 = 40L,
          Var2 = 40L)
)
id_series <- id_series[order(id_series[,1L] * id_series[,2L]),]

mycat("\tRunning sdtw experiments for multiple univariate series\n")
dist_sdtw_univariate_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    benchmarks <- lapply(series, function(this_series) {
        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "sdtw",
                            error.check = FALSE)
            )
        })

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "sdtw_univariate",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# normalized gak univariate
# --------------------------------------------------------------------------------------------------

mycat("\tRunning normalized_gak experiments for multiple univariate series\n")
dist_ngak_univariate_multiple <- dplyr::bind_rows(lapply(num_threads_to_test, function(num_threads) {
    mycat("\t\t")
    RcppParallel::setThreadOptions(num_threads)

    benchmarks <- lapply(series, function(this_series) {
        expressions <- lapply(window_sizes, function(window_size) {
            lapply(1L:nrow(id_series), function(i) {
                bquote(
                    proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                                y = this_series[1L:.(id_series[i, 2L])],
                                method = "gak",
                                window.size = .(window_size),
                                sigma = 100,
                                error.check = FALSE)
                )
            })
        })
        expressions <- unlist(expressions, recursive = FALSE)

        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        cat(".")
        data.frame(distance = "gak_univariate",
                   num_threads = num_threads,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = rep(window_sizes, each = nrow(id_series)),
                   median_time_ms = benchmark$median)
    })

    cat("\n")
    dplyr::bind_rows(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# aggregate
# --------------------------------------------------------------------------------------------------

dist_multiple_results <- dplyr::bind_rows(
    dist_lbk_multiple,
    dist_lbi_multiple,
    dist_dtwlb_multiple,
    dist_sbd_multiple,
    dist_dtw_univariate_multiple,
    dist_dtw_multivariate_multiple,
    dist_sdtw_univariate_multiple,
    dist_ngak_univariate_multiple
)

# Make factor with the given order
dist_multiple_results$distance <- factor(dist_multiple_results$distance,
                                         levels = unique(dist_multiple_results$distance))

# Add some metadata
attr(dist_multiple_results, "proctime") <- proc.time() - t1
attr(dist_multiple_results, "times") <- times

# Clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "dist_multiple_results")))

# ==================================================================================================
# finish
# ==================================================================================================

save("dist_single_results", "dist_multiple_results", file = "dist-results.RData")
cat("\n")
