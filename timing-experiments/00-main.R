suppressPackageStartupMessages({
    library("dtwclust")
    library("microbenchmark")
    library("plyr")
    library("doParallel")
})

#' Set this to TRUE to run a subset of the experiments with less evaluations. "Short" is relative
#' though, it will still take a few hours to complete.
short_experiments <- TRUE

if (short_experiments) message("\nShort experiments active\n") else message("\nShort experiments NOT active\n")

tic <- proc.time()
if (file.exists("read-csv.RData")) load("read-csv.RData") else source("10-read-csv.R")
if (file.exists("dist-results.RData")) load("dist-results.RData") else source("20-distance-experiments.R")
if (file.exists("cent-results.RData")) load("cent-results.RData") else source("30-prototyping-experiments.R")
toc <- proc.time() - tic

message("Finished after: ", toc["elapsed"], " seconds")
