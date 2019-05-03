## ----setup, include = FALSE, cache = FALSE-------------------------------
library("dtwclust")
data("uciCT")

## ----dtw-intuition, out.width = "0.75\\linewidth", fig.cap = "Sample alignment performed by the DTW algorithm between two series. The dashed blue lines exemplify how some points are mapped to each other, which shows how they can be warped in time. Note that the vertical position of each series was artificially altered for visualization."----
dtw_example <- dtw(CharTraj[[1L]], CharTraj[[2L]], keep.internals = TRUE)
plot(dtw_example, type = "two",
     offset = 1, match.indices = 30,
     match.col = "blue",
     xlab = "Time", ylab = "Series")

## ----step-patterns, out.width = "0.45\\linewidth", fig.width = 4, fig.asp = 1, fig.cap = "Two common step patterns used by DTW when traversing the LCM. At each step, the lines denote the allowed directions that can be taken, as well as the weight associated with each one.", fig.subcap = c("\\code{symmetric1} step pattern", "\\code{symmetric2} step pattern")----
plot(symmetric1)
plot(symmetric2)

## ----dtw-window-plot, out.width = "0.6\\linewidth", fig.asp = 1, fig.cap = "Visual representation of the Sakoe-Chiba constraint for DTW. The red elements will not be considered by the algorithm when traversing the LCM."----
dtwWindow.plot(sakoeChibaWindow, window.size = 2, reference = 10, query = 10)

## ----envelope-plot, out.width = "0.75\\linewidth", fig.cap = "Visual representation of a time-series (shown as a solid black line) and its corresponding envelopes based on a Sakoe-Chiba window of size 15. The green dashed line represents the upper envelope, while the red dashed line represents the lower envelope."----
envelopes <- compute_envelope(CharTraj[[2L]], window.size = 15)
matplot(cbind(envelopes$lower, envelopes$upper),
        type = "l", lty = 2, col = 2:3,
        xlab = "Time", ylab = "Series")
lines(CharTraj[[2L]])

## ----using-proxy---------------------------------------------------------
require("TSclust")
proxy::pr_DB$set_entry(FUN = diss.ACF, names = c("ACFD"),
                       loop = TRUE, distance = TRUE,
                       description = "Autocorrelation-based distance")

proxy::dist(CharTraj[3:8], method = "ACFD", upper = TRUE)

## ----benchmarks-against-other-packages-----------------------------------
library(dplyr)
library(TSdist)
library(microbenchmark)
library(ggplot2)

ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(legend.position = "bottom")

#' The dataset has series with different lengths. This threshold specifies the difference there
#' should be between series in the subset that is taken.
length_diff_threshold <- 20L

# Take only two samples of each character
series <- lapply(seq(from = 1L, to = 100L, by = 5L), function(i) { CharTraj[i:(i + 1L)] })
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
series_mv <- lapply(seq(from = 1L, to = 100L, by = 5L), function(i) { CharTrajMV[i:(i + 1L)] })
series_mv <- series_mv[id_ascending]

# Window sizes for the experiments
window_sizes <- seq(from = 10L, to = 50L, by = 10L)
# Number of times each experiment will be repeated (by microbenchmark)
times <- 100L

## ----lbk-single----------------------------------------------------------
lbk_single <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series, function(this_series) {
            # Build expressions to evaluate, substituting window size
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    lb_keogh(x, y, .(window.size), norm = "L2", error.check = FALSE)
                )
            })

            # Extract sample series
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            # Return data frame with results
            data.frame(distance = "LB Keogh (dtwclust)",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        # Bind results for all series and return to global environment
        dplyr::bind_rows(benchmarks)
    }
)

lbk_tsdist_single <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series, function(this_series) {
            # Build expressions to evaluate, substituting window size
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    LBKeoghDistance(y, x, .(window.size * 2L + 1L))
                )
            })

            # Extract sample series
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            # Return data frame with results
            data.frame(distance = "LB Keogh (TSdist)",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        # Bind results for all series and return to global environment
        dplyr::bind_rows(benchmarks)
    }
)

lbk_comparison <- dplyr::bind_rows(lbk_single, lbk_tsdist_single)
lbk_comparison <- dplyr::mutate_if(lbk_comparison, ~ is.character(.) | is.integer(.), as.factor)

ggplot(lbk_comparison,
       aes(x = series_length,
           y = median_time_us,
           group = window_size,
           colour = window_size)) +
    geom_point() +
    geom_line() +
    facet_wrap(~distance, scales = "free_y") +
    scale_color_discrete(name = "Window size") +
    labs(x = "Series' length", y = expression("Median time ("*mu*"s)"))

## ----dtw-single----------------------------------------------------------
dtw_basic_single <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series, function(this_series) {
            # Build expressions to evaluate, substituting window size
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    dtw_basic(x, y, .(window.size), norm = "L2")
                )
            })

            # Extract sample series
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            # Return data frame with results
            data.frame(distance = "dtw_basic (dtwclust)",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        # Bind results for all series and return to global environment
        dplyr::bind_rows(benchmarks)
    }
)

dtw_single <- with(
    new.env(),
    {
        # Loop along extracted subsets
        benchmarks <- lapply(series, function(this_series) {
            # Build expressions to evaluate, substituting window size
            expressions <- lapply(window_sizes, function(window.size) {
                bquote(
                    dtw::dtw(x, y, window.type = "slantedband", window.size = .(window.size), distance.only = TRUE)
                )
            })

            # Extract sample series
            x <- this_series[[1L]]
            y <- this_series[[2L]]

            # Evaluate expressions
            benchmark <- summary(microbenchmark(list = expressions, times = times, unit = "us"))

            # Return data frame with results
            data.frame(distance = "dtw (dtw)",
                       series_length = NROW(x),
                       window_size = window_sizes,
                       median_time_us = benchmark$median,
                       stringsAsFactors = FALSE)
        })

        # Bind results for all series and return to global environment
        dplyr::bind_rows(benchmarks)
    }
)

dtw_comparison <- dplyr::bind_rows(dtw_basic_single, dtw_single)
dtw_comparison <- dplyr::mutate_if(dtw_comparison, ~ is.character(.) | is.integer(.), as.factor)

ggplot(dtw_comparison,
       aes(x = series_length,
           y = median_time_us,
           group = window_size,
           colour = window_size)) +
    geom_point() +
    geom_line() +
    facet_wrap(~distance, scales = "free_y") +
    scale_color_discrete(name = "Window size") +
    labs(x = "Series' length", y = expression("Median time ("*mu*"s)"))

## ----lbk-multiple--------------------------------------------------------
RcppParallel::setThreadOptions(1L)
series <- reinterpolate(CharTraj, 100L)
window_size <- 30L
id_series <- cbind(c(seq(from = 10L, to = 40L, by = 10L),
                     rep(50L, 6L)),
                   seq(from = 10L, to = 100L, by = 10L))

lbk_multiple <- {
    # Build expressions for proxy::dist
    expressions <- lapply(1L:nrow(id_series), function(i) {
        bquote(
            proxy::dist(x = series[1L:.(id_series[i, 1L])],
                        y = series[1L:.(id_series[i, 2L])],
                        method = "lbk",
                        norm = "L2",
                        window.size = .(window_size),
                        error.check = FALSE)
        )
    })

    # Evaluate expressions
    benchmark <- summary(microbenchmark(list = expressions, times = 30L, unit = "ms"))

    # Return data frame with results
    data.frame(distance = "LB Keogh (dtwclust)",
               num_x = id_series[,1L],
               num_y = id_series[,2L],
               num_total = id_series[,1L] * id_series[,2L],
               series_length = NROW(series[[1L]]),
               window_size = window_size,
               median_time_ms = benchmark$median,
               stringsAsFactors = FALSE)
}

lbk_tsdist_multiple <- {
    # Build expressions for proxy::dist
    expressions <- lapply(1L:nrow(id_series), function(i) {
        # x and y inverted because the TSdist version calculates envelops for x
        bquote(
            proxy::dist(y = series[1L:.(id_series[i, 1L])],
                        x = series[1L:.(id_series[i, 2L])],
                        method = "TSDistances",
                        distance = "lb.keogh",
                        window.size = .(window_size * 2L + 1L))
        )
    })

    # Evaluate expressions
    benchmark <- summary(microbenchmark(list = expressions, times = 3L, unit = "ms"))

    # Return data frame with results
    data.frame(distance = "LB Keogh (TSdist)",
               num_x = id_series[,1L],
               num_y = id_series[,2L],
               num_total = id_series[,1L] * id_series[,2L],
               series_length = NROW(series[[1L]]),
               window_size = window_size,
               median_time_ms = benchmark$median,
               stringsAsFactors = FALSE)
}

lbk_comparison <- dplyr::bind_rows(lbk_multiple, lbk_tsdist_multiple)
lbk_comparison$vbreaks <- 2500L

ggplot(lbk_comparison,
       aes(x = num_total,
           y = median_time_ms,
           colour = num_y)) +
    geom_point(size = 3) +
    geom_line() +
    geom_vline(aes(xintercept = vbreaks), colour = "black", linetype = "longdash") +
    facet_wrap(~ distance, scales = "free_y") +
    scale_color_continuous(name = "Amount of warping envelopes needed") +
    labs(x = "Total number of distance calculations", y = "Median time (ms)")

## ----dtw-multiple--------------------------------------------------------
series <- CharTraj
id_series <- matrix(rep(seq(from = 10L, to = 50L, by = 10L), 2L), ncol = 2L)

dtw_basic_multiple <- {
    # Build expressions for proxy::dist
    expressions <- lapply(1L:nrow(id_series), function(i) {
        bquote(
            proxy::dist(x = series[1L:.(id_series[i, 1L])],
                        y = series[1L:.(id_series[i, 2L])],
                        method = "dtw_basic",
                        norm = "L2",
                        window.size = .(window_size))
        )
    })

    # Evaluate expressions
    benchmark <- summary(microbenchmark(list = expressions, times = 10L, unit = "s"))

    # Return data frame with results
    data.frame(distance = "dtw_basic (dtwclust)",
               num_total = id_series[,1L] * id_series[,2L],
               series_length = NROW(series[[1L]]),
               window_size = window_size,
               median_time_s = benchmark$median,
               stringsAsFactors = FALSE)
}

dtw_multiple <- {
    # Build expressions for proxy::dist
    expressions <- lapply(1L:nrow(id_series), function(i) {
        bquote(
            proxy::dist(y = series[1L:.(id_series[i, 1L])],
                        x = series[1L:.(id_series[i, 2L])],
                        method = "dtw",
                        window.type = "slantedband",
                        window.size = .(window_size))
        )
    })

    # Evaluate expressions
    benchmark <- summary(microbenchmark(list = expressions, times = 10L, unit = "s"))

    # Return data frame with results
    data.frame(distance = "dtw (dtw)",
               num_total = id_series[,1L] * id_series[,2L],
               series_length = NROW(series[[1L]]),
               window_size = window_size,
               median_time_s = benchmark$median,
               stringsAsFactors = FALSE)
}

dtw_comparison <- dplyr::bind_rows(dtw_basic_multiple, dtw_multiple)

ggplot(dtw_comparison,
       aes(x = num_total,
           y = median_time_s)) +
    geom_point(size = 3) +
    geom_line() +
    facet_wrap(~ distance, scales = "free_y") +
    labs(x = "Total number of distance calculations", y = "Median time (s)")

RcppParallel::setThreadOptions()

## ----sbd-alignment, out.width = "0.45\\linewidth", fig.width = 6, fig.cap = "Visualization of the NCCc-based alignment performed on two sample series. After alignment, the second (red) series is either truncated and/or prepended/appended with zeros so that its length matches the first(black) series.", fig.subcap = c("Series before alignment", "Series after alignment")----
matplot(cbind(CharTraj[[61L]], CharTraj[[65L]]),
        type = "l", lty = 1L,
        xlab = "Time", ylab = "Series")

sbd_align <- SBD(CharTraj[[61L]], CharTraj[[65L]])

matplot(cbind(CharTraj[[61L]], sbd_align$yshift),
        type = "l", lty = 1L,
        xlab = "Time", ylab = "Series")

## ----example-pc, echo = TRUE, warning = FALSE----------------------------
# Linear reinterpolation to same length
data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
# z-normalization
data <- zscore(data[60L:100L])

pc_dtw <- tsclust(data, k = 4L, seed = 8L,
                  distance = "dtw_basic", centroid = "dba",
                  norm = "L2", window.size = 20L)

pc_ks <- tsclust(data, k = 4L, seed = 8L,
                 distance = "sbd", centroid = "shape")

pc_tp <- tsclust(data, k = 4L, type = "tadpole", seed = 8L,
                 control = tadpole_control(dc = 1.5, window.size = 20L))

sapply(list(DTW = pc_dtw, kShape = pc_ks, TADPole = pc_tp),
       cvi, b = CharTrajLabels[60L:100L], type = "VI")

## ----example-custom-cent, echo = TRUE------------------------------------
weighted_mean_cent <- function(x, cl_id, k, cent, cl_old, ..., weights) {
    x <- Map(x, weights, f = function(ts, w) { w * ts })
    x_split <- split(x, cl_id)
    new_cent <- lapply(x_split, function(xx) {
        xx <- do.call(rbind, xx)
        colMeans(xx)
    })
}

data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
weights <- rep(c(0.9,1.1), each = 5L)
tsclust(data[1L:10L], type = "p", k = 2L,
        distance = "Manhattan",
        centroid = weighted_mean_cent,
        seed = 123,
        args = tsclust_args(cent = list(weights = weights)))

## ----example-cvis, echo = TRUE-------------------------------------------
data <- CharTraj[1L:20L]
pc_k <- tsclust(data, k = 3L:5L,
                distance = "dtw_basic", centroid = "pam",
                seed = 94L)
names(pc_k) <- paste0("k_", 3L:5L)
sapply(pc_k, cvi, type = "internal")

## ----example-clue-p, echo = TRUE-----------------------------------------
require("clue")

pc_4 <- tsclust(data, type = "p", k = 4L,
                distance = "dtw_basic", centroid = "pam",
                control = partitional_control(nrep = 5L),
                seed = 95L)

names(pc_4) <- paste0("r_", 1L:5L)
pc_4 <- cl_ensemble(list = pc_4)
cl_dissimilarity(pc_4)

# Confusion matrix
table(Medoid = cl_class_ids(cl_medoid(pc_4)),
      "True Classes" = rep(c(4L, 3L, 1L, 2L), each = 5L))

## ----compare-clusterings-------------------------------------------------
require("doParallel")
workers <- makeCluster(detectCores())
invisible(clusterEvalQ(workers, library(dtwclust)))
registerDoParallel(workers)

cfg <- compare_clusterings_configs(
    types = "partitional",
    k = 20L,
    controls = list(
        partitional = partitional_control(
            iter.max = 20L
        )
    ),
    distances = pdc_configs(
        "distance",
        partitional = list(
            dtw_basic = list(
                window.size = seq(from = 10L, to = 30L, by = 5L),
                norm = c("L1", "L2")
            )
        )
    ),
    centroids = pdc_configs(
        "centroid",
        share.config = c("p"),
        dba = list(
            window.size = seq(from = 10L, to = 30L, by = 5L),
            norm = c("L1", "L2")
        )
    ),
    no.expand = c(
        "window.size",
        "norm"
    )
)

evaluators <- cvi_evaluators("ARI", ground.truth = CharTrajLabels)

comparison <- compare_clusterings(CharTraj, types = "partitional",
                                  configs = cfg, seed = 8L,
                                  score.clus = evaluators$score,
                                  pick.clus = evaluators$pick)

stopCluster(workers); registerDoSEQ()

# some rows and columns from the results data frame
head(comparison$results$partitional[, c("config_id", "distance", "centroid",
                                        "window.size_distance", "norm_distance",
                                        "ARI")])

clusters <- repeat_clustering(CharTraj, comparison, comparison$pick$config_id)

matrix(clusters@cluster, ncol = 5L, byrow = TRUE)
