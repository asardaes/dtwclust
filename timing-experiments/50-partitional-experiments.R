# ==================================================================================================
# TADPole experiments
# ==================================================================================================

existing_objects <- ls(all.names = TRUE)

# --------------------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------------------

# Take only one variable of each character
series <- univariate_series[seq(from = 1L, to = length(univariate_series), by = 3L)]
# Length of reinterpolated series
new_length <- median(sapply(series, function(s) { length(s[[1L]]) }))

if (short_experiments) {
    # How many of **each** character to consider
    num_series <- seq(from = 10L, to = 30L, by = 10L)
    # Window sizes to consider
    window_sizes <- seq(from = 40L, to = 60L, by = 20L)
    # Number of repetitions for the PAM vs repetitions experiment
    repetitions <- c(1L, 5L, 10L)
    # Number of clusters to test for sparse PAM case
    sparse_k <- seq(from = 5L, to = 20L, by = 5L)
    # Number of evaluations
    times <- 2L

} else {
    # How many of **each** character to consider
    num_series <- seq(from = 5L, to = 30L, by = 5L)
    # Window sizes to consider
    window_sizes <- seq(from = 20L, to = 100L, by = 20L)
    # Number of repetitions for the PAM vs repetitions experiment
    repetitions <- 1L:10L
    # Number of clusters to test for sparse PAM case
    sparse_k <- seq(from = 2L, to = 18L, by = 4L)
    # Number of evaluations
    times <- 15L
}

# NOTE: all clustering experiments will use tsclust() to include overhead of corresponding checks

t1 <- proc.time()

# ==================================================================================================
# dtw_basic vs dtw_lb
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# PAM
# --------------------------------------------------------------------------------------------------

cat("\tRunning dtw_basic vs dtw_lb clustering experiments (PAM)\n")
clus_dtwb_dtwlb_pam_results <- dplyr::bind_rows(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset and reinterpolate to equal length
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)
    series <- reinterpolate(series, new.length = new_length)

    benchmarks <- dplyr::bind_rows(lapply(window_sizes, function(window_size) {
        times <- sapply(1L:times, function(dummy) {
            tsc_dtwb <- tsclust(series = series, k = 20L, type = "partitional",
                                distance = "dtw_basic", centroid = "pam",
                                seed = window_size, trace = FALSE, error.check = FALSE,
                                control = partitional_control(pam.precompute = TRUE,
                                                              iter.max = 10L),
                                args = tsclust_args(dist = list(window.size = window_size,
                                                                norm = "L1",
                                                                step.pattern = dtw::symmetric1)))

            tsc_dtwbs <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_basic", centroid = "pam",
                                 seed = window_size, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(pam.precompute = FALSE,
                                                               pam.sparse = TRUE,
                                                               iter.max = 10L),
                                 args = tsclust_args(dist = list(window.size = window_size,
                                                                 norm = "L1",
                                                                 step.pattern = dtw::symmetric1)))

            tsc_dtwlb <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_lb", centroid = "pam",
                                 seed = window_size, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(pam.precompute = FALSE,
                                                               pam.sparse = FALSE,
                                                               iter.max = 10L),
                                 args = tsclust_args(dist = list(window.size = window_size,
                                                                 norm = "L1",
                                                                 step.pattern = dtw::symmetric1)))

            c(dtwb = tsc_dtwb@proctime[["elapsed"]],
              dtwbs = tsc_dtwbs@proctime[["elapsed"]],
              dtwlb = tsc_dtwlb@proctime[["elapsed"]])
        })

        median_times <- apply(times, 1L, median)

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = 20L,
                   window_size = window_size,
                   dtw_basic_median_time_s = median_times[["dtwb"]],
                   dtw_basic_sparse_median_time_s = median_times[["dtwbs"]],
                   dtw_lb_median_time_s = median_times[["dtwlb"]])
    }))

    cat("\n")
    benchmarks
}))

# --------------------------------------------------------------------------------------------------
# PAM vs repetitions
# --------------------------------------------------------------------------------------------------

cat("\tRunning dtw_basic vs dtw_lb clustering experiments (PAM vs nrep)\n")
clus_dtwb_dtwlb_pamrep_results <- dplyr::bind_rows(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset and reinterpolate to equal length
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)
    series <- reinterpolate(series, new.length = new_length)

    benchmarks <- dplyr::bind_rows(lapply(repetitions, function(nrep) {
        times <- sapply(1L:times, function(dummy) {
            tsc_dtwb <- tsclust(series = series, k = 20L, type = "partitional",
                                distance = "dtw_basic", centroid = "pam",
                                seed = nrep, trace = FALSE, error.check = FALSE,
                                control = partitional_control(pam.precompute = TRUE,
                                                              iter.max = 10L,
                                                              nrep = nrep),
                                args = tsclust_args(dist = list(window.size = 20L,
                                                                norm = "L1",
                                                                step.pattern = dtw::symmetric1)))

            tsc_dtwbs <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_basic", centroid = "pam",
                                 seed = nrep, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(pam.precompute = FALSE,
                                                               pam.sparse = TRUE,
                                                               iter.max = 10L,
                                                               nrep = nrep),
                                 args = tsclust_args(dist = list(window.size = 20L,
                                                                 norm = "L1",
                                                                 step.pattern = dtw::symmetric1)))

            tsc_dtwlb <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_lb", centroid = "pam",
                                 seed = nrep, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(pam.precompute = FALSE,
                                                               pam.sparse = FALSE,
                                                               iter.max = 10L,
                                                               nrep = nrep),
                                 args = tsclust_args(dist = list(window.size = 20L,
                                                                 norm = "L1",
                                                                 step.pattern = dtw::symmetric1)))

            distmat <- if (nrep > 1L) tsc_dtwbs[[1L]]@distmat else tsc_dtwbs@distmat

            c(dtwb = if (nrep > 1L) tsc_dtwb[[1L]]@proctime[["elapsed"]] else tsc_dtwb@proctime[["elapsed"]],
              dtwbs = if (nrep > 1L) tsc_dtwbs[[1L]]@proctime[["elapsed"]] else tsc_dtwbs@proctime[["elapsed"]],
              dtwlb = if (nrep > 1L) tsc_dtwlb[[1L]]@proctime[["elapsed"]] else tsc_dtwlb@proctime[["elapsed"]],
              sparse_distmat_filled = 100 * sum(distmat != 0) / length(distmat))
        })

        median_times <- apply(times, 1L, median) # sparse_distmat_filled should not vary

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = 20L,
                   num_repetitions = nrep,
                   dtw_basic_median_time_s = median_times[["dtwb"]],
                   dtw_basic_sparse_median_time_s = median_times[["dtwbs"]],
                   dtw_lb_median_time_s = median_times[["dtwlb"]],
                   sparse_distmat_filled_percent = median_times[["sparse_distmat_filled"]])
    }))

    cat("\n")
    benchmarks
}))

# --------------------------------------------------------------------------------------------------
# DBA
# --------------------------------------------------------------------------------------------------

cat("\tRunning dtw_basic vs dtw_lb clustering experiments (DBA)\n")
clus_dtwb_dtwlb_dba_results <- dplyr::bind_rows(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset and reinterpolate to equal length
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)
    series <- reinterpolate(series, new.length = new_length)

    benchmarks <- dplyr::bind_rows(lapply(window_sizes, function(window_size) {
        times <- sapply(1L:times, function(dummy) {
            tsc_dtwb <- tsclust(series = series, k = 20L, type = "partitional",
                                distance = "dtw_basic", centroid = "dba",
                                seed = window_size, trace = FALSE, error.check = FALSE,
                                control = partitional_control(iter.max = 10L),
                                args = tsclust_args(dist = list(window.size = window_size,
                                                                norm = "L1",
                                                                step.pattern = dtw::symmetric1),
                                                    cent = list(window.size = window_size,
                                                                max.iter = 15L,
                                                                step.pattern = dtw::symmetric1)))

            tsc_dtwlb <- tsclust(series = series, k = 20L, type = "partitional",
                                 distance = "dtw_lb", centroid = "dba",
                                 seed = window_size, trace = FALSE, error.check = FALSE,
                                 control = partitional_control(iter.max = 10L),
                                 args = tsclust_args(dist = list(window.size = window_size,
                                                                 norm = "L1",
                                                                 step.pattern = dtw::symmetric1),
                                                     cent = list(window.size = window_size,
                                                                 max.iter = 15L,
                                                                 step.pattern = dtw::symmetric1)))

            c(dtwb = tsc_dtwb@proctime[["elapsed"]],
              dtwlb = tsc_dtwlb@proctime[["elapsed"]])
        })

        median_times <- apply(times, 1L, median)

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = 20L,
                   window_size = window_size,
                   dtw_basic_median_time_s = median_times[["dtwb"]],
                   dtw_lb_median_time_s = median_times[["dtwlb"]])
    }))

    cat("\n")
    benchmarks
}))

# ==================================================================================================
# sparse PAM vs different k values
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# non-symmetric
# --------------------------------------------------------------------------------------------------

cat("\tRunning experiments for sparse PAM vs different k \n")
clus_pam_sparse_k_results <- dplyr::bind_rows(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)

    benchmarks <- dplyr::bind_rows(lapply(sparse_k, function(k) {
        times <- sapply(1L:times, function(dummy) {
            tsc_nonsparse <- tsclust(series = series, k = k, type = "partitional",
                                     distance = "dtw_basic", centroid = "pam",
                                     seed = k, trace = FALSE, error.check = FALSE,
                                     control = partitional_control(pam.precompute = TRUE,
                                                                   pam.sparse = FALSE,
                                                                   iter.max = 10L),
                                     args = tsclust_args(dist = list(window.size = 20L,
                                                                     norm = "L1",
                                                                     step.pattern = dtw::symmetric1)))

            tsc_sparse <- tsclust(series = series, k = k, type = "partitional",
                                  distance = "dtw_basic", centroid = "pam",
                                  seed = k, trace = FALSE, error.check = FALSE,
                                  control = partitional_control(pam.precompute = FALSE,
                                                                pam.sparse = TRUE,
                                                                iter.max = 10L),
                                  args = tsclust_args(dist = list(window.size = 20L,
                                                                  norm = "L1",
                                                                  step.pattern = dtw::symmetric1)))

            distmat <- tsc_sparse@distmat

            c(non_sparse = tsc_nonsparse@proctime[["elapsed"]],
              sparse = tsc_sparse@proctime[["elapsed"]],
              sparse_distmat_filled = 100 * sum(distmat != 0) / length(distmat))
        })

        median_times <- apply(times, 1L, median) # sparse_distmat_filled should not vary

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = k,
                   non_sparse_median_time_s = median_times[["non_sparse"]],
                   sparse_median_time_s = median_times[["sparse"]],
                   sparse_distmat_filled_percent = median_times[["sparse_distmat_filled"]])
    }))

    cat("\n")
    benchmarks
}))

# --------------------------------------------------------------------------------------------------
# symmetric
# --------------------------------------------------------------------------------------------------

cat("\tRunning experiments for sparse, symmetric PAM vs different k \n")
clus_pam_sparse_symmetric_k_results <- dplyr::bind_rows(lapply(num_series, function(num_series) {
    cat("\t\t")

    # Get subset
    series <- lapply(series, function(s) { s[1L:num_series] })
    series <- unlist(series, recursive = FALSE)

    benchmarks <- dplyr::bind_rows(lapply(sparse_k, function(k) {
        times <- sapply(1L:times, function(dummy) {
            tsc_nonsparse <- tsclust(series = series, k = k, type = "partitional",
                                     distance = "sbd", centroid = "pam",
                                     seed = k, trace = FALSE, error.check = FALSE,
                                     control = partitional_control(pam.precompute = TRUE,
                                                                   pam.sparse = FALSE,
                                                                   iter.max = 10L))

            tsc_sparse <- tsclust(series = series, k = k, type = "partitional",
                                  distance = "sbd", centroid = "pam",
                                  seed = k, trace = FALSE, error.check = FALSE,
                                  control = partitional_control(pam.precompute = FALSE,
                                                                pam.sparse = TRUE,
                                                                iter.max = 10L))

            distmat <- tsc_sparse@distmat

            c(non_sparse = tsc_nonsparse@proctime[["elapsed"]],
              sparse = tsc_sparse@proctime[["elapsed"]],
              sparse_distmat_filled = 100 * sum(distmat != 0) / length(distmat))
        })

        median_values <- apply(times, 1L, median)

        cat(".")
        # Return data frame
        data.frame(num_series = length(series),
                   k = k,
                   non_sparse_median_time_s = median_values[["non_sparse"]],
                   sparse_median_time_s = median_values[["sparse"]],
                   sparse_distmat_filled_percent = median_values[["sparse_distmat_filled"]])
    }))

    cat("\n")
    benchmarks
}))

# ==================================================================================================
# aggregate
# ==================================================================================================

partitional_results <- list(
    dtwlb_vs_dtwbasic = list(
        pam = clus_dtwb_dtwlb_pam_results,
        pam_vs_reps = clus_dtwb_dtwlb_pamrep_results,
        dba = clus_dtwb_dtwlb_dba_results
    ),
    sparse_pam_k = list(
        non_symmetric = clus_pam_sparse_k_results,
        symmetric = clus_pam_sparse_symmetric_k_results
    )
)

# Add some metadata
attr(partitional_results, "proctime") <- proc.time() - t1
attr(partitional_results, "times") <- times

# ==================================================================================================
# finish
# ==================================================================================================

# Clean
rm(list = setdiff(ls(all.names = TRUE), c(existing_objects, "partitional_results")))
save("partitional_results", file = "partitional-results.RData")
cat("\n")
