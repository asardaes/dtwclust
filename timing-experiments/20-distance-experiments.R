# ==================================================================================================
# get sample series and sort them by length
# ==================================================================================================

existing_objects <- ls(all.names = TRUE)

series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
series <- lapply(series, function(s) { s[1L:2L] })
len <- sapply(series, function(s) { lengths(s)[1L] })
id_ascending <- sort(len, decreasing = FALSE, index.return = TRUE)
series <- series[id_ascending$ix][!duplicated(id_ascending$x)]
normalized_series <- zscore(series)
series_mv <- lapply(multivariate_series, function(s) { s[1L:2L] })
series_mv <- series_mv[id_ascending$ix][!duplicated(id_ascending$x)]

# ==================================================================================================
# distance experiments using single series (1 to 1)
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# lb_keogh
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "lb_keogh",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# lb_improved
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "lb_improved",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# sbd
# --------------------------------------------------------------------------------------------------

dist_sbd_single <- with(
    new.env(),
    {
        benchmarks <- lapply(normalized_series, function(this_series) {
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            benchmark <- summary(microbenchmark(SBD(x, y, error.check = FALSE, return.shifted = FALSE),
                                                times = 100L, unit = "us"))

            data.frame(distance = "sbd",
                       series_length = NROW(x),
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# dtw univariate
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "dtw_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# dtw_multivariate
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "dtw_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# unnormalized gak univariate
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "unnormalized_gak_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# unnormalized gak multivariate
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "unnormalized_gak_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# normalized gak univariate
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "normalized_gak_univariate",
                       series_length = NROW(x),
                       window_size = c(window_sizes, NA),
                       timing_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        plyr::rbind.fill(benchmarks)
    }
)

# --------------------------------------------------------------------------------------------------
# normalized gak multivariate
# --------------------------------------------------------------------------------------------------

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

            benchmark <- summary(microbenchmark(list = expressions, times = 100L, unit = "us"))

            data.frame(distance = "normalized_gak_multivariate",
                       series_length = nrow(x),
                       window_size = c(window_sizes, NA),
                       timing_us = benchmark$median,
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

# ==================================================================================================
# finish
# ==================================================================================================

## temporary
ggplot(dist_single_results,
       aes(x = factor(series_length),
           y = timing_us,
           group = factor(window_size),
           colour = factor(window_size))) +
    geom_line() +
    facet_wrap(~distance, scales = "free_y") +
    theme_bw()

## clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "dist_single_results")))
save("dist_single_results", file = "dist-results.RData")
