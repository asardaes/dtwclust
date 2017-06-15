library("dtwclust")
library("microbenchmark")
library("plyr")
library("doParallel")

short_experiments <- TRUE
length_diff_threshold <- 20L

if (file.exists("read-csv.RData")) load("read-csv.RData") else source("10-read-csv.R")
if (file.exists("dist-results.RData")) load("dist-results.RData") else source("20-distance-experiments.R")
