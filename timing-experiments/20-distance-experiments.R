existing_objects <- ls(all.names = TRUE)

# ==================================================================================================
# distance experiments using single series (1 to 1)
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# get sample series and sort them by length
# --------------------------------------------------------------------------------------------------

series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
series <- lapply(series, function(s) { s[1L:2L] })
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
series_mv <- lapply(multivariate_series, function(s) { s[1L:2L] })
series_mv <- series_mv[id_ascending]

# --------------------------------------------------------------------------------------------------
# lb_keogh
# --------------------------------------------------------------------------------------------------

cat("\tRunning lb_keogh experiments for single series\n")
dist_lbk_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    lb_keogh(x, y, .(window.size), error.check = FALSE)
                )
            })

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "lb_keogh",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# lb_improved
# --------------------------------------------------------------------------------------------------

cat("\tRunning lb_improved experiments for single series\n")
dist_lbi_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    lb_improved(x, y, .(window.size), error.check = FALSE)
                )
            })

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "lb_improved",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# sbd
# --------------------------------------------------------------------------------------------------

cat("\tRunning sbd experiments for single series\n")
dist_sbd_single <- with(
    new.env(),
    {
        benchmarks <- lapply(series, function(this_series) {
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(
                SBD(x, y, error.check = FALSE, return.shifted = FALSE),
                times = times, unit = "us"
            ))

            data.frame(distance = "sbd",
                       series_length = NROW(x),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# dtw univariate
# --------------------------------------------------------------------------------------------------

cat("\tRunning dtw experiments for single univariate series\n")
dist_dtw_univariate_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

        benchmarks <- lapply(series, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    dtw_basic(x, y, .(window.size), error.check = FALSE)
                )
            })

            expressions <- c(expressions, list(bquote(dtw_basic(x, y, error.check = FALSE))))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "dtw_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# dtw_multivariate
# --------------------------------------------------------------------------------------------------

cat("\tRunning dtw experiments for single multivariate series\n")
dist_dtw_multivariate_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

        benchmarks <- lapply(series_mv, function(this_series) {
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    dtw_basic(x, y, .(window.size), error.check = FALSE)
                )
            })

            expressions <- c(expressions, list(bquote(dtw_basic(x, y, error.check = FALSE))))

            x <- this_series[[1L]]
            y <- this_series[[2L]]

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "dtw_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# unnormalized gak univariate
# --------------------------------------------------------------------------------------------------

cat("\tRunning unnormalized_gak experiments for single univariate series\n")
dist_unnormalized_gak_univariate_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

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

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "unnormalized_gak_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# unnormalized gak multivariate
# --------------------------------------------------------------------------------------------------

cat("\tRunning unnormalized_gak experiments for single multivariate series\n")
dist_unnormalized_gak_multivariate_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

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

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "unnormalized_gak_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# normalized gak univariate
# --------------------------------------------------------------------------------------------------

cat("\tRunning normalized_gak experiments for single univariate series\n")
dist_normalized_gak_univariate_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

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

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "normalized_gak_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# normalized gak multivariate
# --------------------------------------------------------------------------------------------------

cat("\tRunning normalized_gak experiments for single multivariate series\n")
dist_normalized_gak_multivariate_single <- with(
    new.env(),
    {
        window_sizes <- seq(from = 10L, to = 100L, by = 10L)

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

            times <- if (short_experiments) 10L else 100L
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            data.frame(distance = "normalized_gak_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# aggregate
# --------------------------------------------------------------------------------------------------

dist_single_results <- plyr::rbind.fill(
    dist_lbk_single,
    dist_lbi_single,
    dist_sbd_single,
    dist_dtw_univariate_single,
    dist_dtw_multivariate_single,
    dist_unnormalized_gak_univariate_single,
    dist_unnormalized_gak_multivariate_single,
    dist_normalized_gak_univariate_single,
    dist_normalized_gak_multivariate_single
)

## make factor with the given order
dist_single_results$distance <- factor(dist_single_results$distance,
                                       levels = unique(dist_single_results$distance))

## clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "dist_single_results")))
existing_objects <- ls(all.names = TRUE)

# ==================================================================================================
# distance experiments using multiple series
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# get series and sort them by length
# --------------------------------------------------------------------------------------------------

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
num_workers_to_test <- 1L:4L

if (short_experiments) {
    series <- series[1L:2L]
    series_mv <- series_mv[1L:2L]
    num_workers_to_test <- c(2L, 4L)
}

cat("\n")

# --------------------------------------------------------------------------------------------------
# lb_keogh
# --------------------------------------------------------------------------------------------------

cat("\tRunning lb_keogh experiments for multiple series\n")
dist_lbk_multiple <- plyr::rbind.fill(lapply(num_workers_to_test, function(num_workers) {
    registerDoParallel(workers <- makeCluster(num_workers))
    invisible(clusterEvalQ(workers, library("dtwclust")))
    invisible(clusterEvalQ(workers, library("microbenchmark")))

    window_size <- 50L
    benchmarks <- lapply(series, function(this_series) {
        if (short_experiments) {
            id_series <- cbind(seq(from = 10L, to = 100L, by = 10L),
                               seq(from = 10L, to = 100L, by = 10L))
        } else {
            id_series <- rbind(
                expand.grid(seq(from = 10L, to = 100L, by = 10L), 10L),
                expand.grid(100L, seq(from = 20L, to = 100L, by = 10L)),
                cbind(Var1 = seq(from = 20L, to = 90L, by = 10L),
                      Var2 = seq(from = 20L, to = 90L, by = 10L))
            )

            id_series <- id_series[order(id_series[,1L] * id_series[,2L]),]
        }

        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "lbk",
                            window.size = .(window_size),
                            error.check = FALSE)
            )
        })

        times <- if (short_experiments) 10L else 50L
        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        data.frame(distance = "lb_keogh",
                   num_workers = num_workers,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = window_size,
                   median_time_ms = benchmark$median,
                   stringsAsFactors = FALSE)
    })

    stopCluster(workers)
    registerDoSEQ()
    rm(workers)

    plyr::rbind.fill(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# lb_improved
# --------------------------------------------------------------------------------------------------

cat("\tRunning lb_improved experiments for multiple series\n")
dist_lbi_multiple <- plyr::rbind.fill(lapply(num_workers_to_test, function(num_workers) {
    registerDoParallel(workers <- makeCluster(num_workers))
    invisible(clusterEvalQ(workers, library("dtwclust")))
    invisible(clusterEvalQ(workers, library("microbenchmark")))

    window_size <- 50L
    benchmarks <- lapply(series, function(this_series) {
        if (short_experiments) {
            id_series <- cbind(seq(from = 10L, to = 100L, by = 10L),
                               seq(from = 10L, to = 100L, by = 10L))
        } else {
            id_series <- rbind(
                expand.grid(seq(from = 10L, to = 100L, by = 10L), 10L),
                expand.grid(100L, seq(from = 20L, to = 100L, by = 10L)),
                cbind(Var1 = seq(from = 20L, to = 90L, by = 10L),
                      Var2 = seq(from = 20L, to = 90L, by = 10L))
            )

            id_series <- id_series[order(id_series[,1L] * id_series[,2L]),]
        }

        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "lbi",
                            window.size = .(window_size),
                            error.check = FALSE)
            )
        })

        times <- if (short_experiments) 10L else 50L
        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        data.frame(distance = "lb_improved",
                   num_workers = num_workers,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   window_size = window_size,
                   median_time_ms = benchmark$median,
                   stringsAsFactors = FALSE)
    })

    stopCluster(workers)
    registerDoSEQ()
    rm(workers)

    plyr::rbind.fill(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# sbd
# --------------------------------------------------------------------------------------------------

cat("\tRunning sbd experiments for multiple series\n")
dist_sbd_multiple <- plyr::rbind.fill(lapply(num_workers_to_test, function(num_workers) {
    registerDoParallel(workers <- makeCluster(num_workers))
    invisible(clusterEvalQ(workers, library("dtwclust")))
    invisible(clusterEvalQ(workers, library("microbenchmark")))

    benchmarks <- lapply(series, function(this_series) {
        if (short_experiments) {
            id_series <- cbind(seq(from = 10L, to = 100L, by = 10L),
                               seq(from = 10L, to = 100L, by = 10L))
        } else {
            id_series <- rbind(
                expand.grid(seq(from = 10L, to = 100L, by = 10L), 10L),
                expand.grid(100L, seq(from = 20L, to = 100L, by = 10L)),
                cbind(Var1 = seq(from = 20L, to = 90L, by = 10L),
                      Var2 = seq(from = 20L, to = 90L, by = 10L))
            )

            id_series <- id_series[order(id_series[,1L] * id_series[,2L]),]
        }

        expressions <- lapply(1L:nrow(id_series), function(i) {
            bquote(
                proxy::dist(x = this_series[1L:.(id_series[i, 1L])],
                            y = this_series[1L:.(id_series[i, 2L])],
                            method = "sbd",
                            error.check = FALSE)
            )
        })

        times <- if (short_experiments) 10L else 50L
        benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "ms"))

        data.frame(distance = "sbd",
                   num_workers = num_workers,
                   num_x = id_series[,1L],
                   num_y = id_series[,2L],
                   num_total = id_series[,1L] * id_series[,2L],
                   series_length = NROW(this_series[[1L]]),
                   median_time_ms = benchmark$median,
                   stringsAsFactors = FALSE)
    })

    stopCluster(workers)
    registerDoSEQ()
    rm(workers)

    plyr::rbind.fill(benchmarks)
}))

# --------------------------------------------------------------------------------------------------
# aggregate
# --------------------------------------------------------------------------------------------------

dist_multiple_results <- plyr::rbind.fill(
    dist_lbk_multiple,
    dist_lbi_multiple,
    dist_sbd_multiple
)

## make factor with the given order
dist_multiple_results$distance <- factor(dist_multiple_results$distance,
                                         levels = unique(dist_multiple_results$distance))

## clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "dist_multiple_results")))

# ==================================================================================================
# finish
# ==================================================================================================

## temporary
ggplot(dist_single_results,
       aes(x = factor(series_length),
           y = median_time_us,
           group = factor(window_size),
           colour = factor(window_size))) +
    geom_line() +
    facet_wrap(~distance, scales = "free_y") +
    theme_bw()

ggplot(dist_multiple_results,
       aes(x = num_total,
           y = median_time_ms,
           colour = num_y,
           shape = factor(series_length))) +
    geom_point(size = 3) +
    facet_grid(distance ~ num_workers) +
    theme_bw()

# save("dist_single_results", "dist_multiple_results", file = "dist-results.RData")
