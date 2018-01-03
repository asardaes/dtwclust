suppressPackageStartupMessages({
    library("dtwclust")
    library("microbenchmark")
    library("dplyr")
    library("doParallel")
})

#' Set this to TRUE to run a subset of the experiments with less evaluations. "Short" is relative
#' though, it will still take a few hours to complete.
#' The short experiments were used during initial setup to fine-tune the parameters.
short_experiments <- FALSE

if (short_experiments) message("\nShort experiments active\n") else message("\nShort experiments NOT active\n")

tic <- proc.time()
if (file.exists("read-csv.RData")) load("read-csv.RData") else source("10-read-csv.R")
if (file.exists("dist-results.RData")) load("dist-results.RData") else source("20-distance-experiments.R")
if (file.exists("cent-results.RData")) load("cent-results.RData") else source("30-prototyping-experiments.R")
if (file.exists("tadpole-results.RData")) load("tadpole-results.RData") else source("40-tadpole-experiments.R")
if (file.exists("partitional-results.RData")) load("partitional-results.RData") else source("50-partitional-experiments.R")
toc <- proc.time() - tic

dtwclustTimings <- list(
    dist = list(
        single = dist_single_results,
        multiple = dist_multiple_results
    ),
    cent = cent_results,
    tadpole = clus_tadpole_results,
    partitional = partitional_results
)

file <- if (short_experiments) "dtwclustTimings.RData" else "../data/dtwclustTimings.rda"
save("dtwclustTimings", file = file)

message("\nFinished after: ", toc["elapsed"], " seconds")
