#### Load data
data(uciCT)

# Reinterpolate to same length and coerce as matrix
data <- t(sapply(CharTraj, reinterpolate, newLength = 180))

# Subset for speed
data <- data[1:20, ]
labels <- CharTrajLabels[1:20]

# Controls
ctrl <- list(trace = TRUE, window.size = 20L)

#### Simple partitional clustering with L2 distance and PAM
kc.l2 <- dtwclust(data, k = 4, distance = "L2", centroid = "pam",
                  seed = 3247, control = ctrl)
cat("Rand index for L2+PAM:", randIndex(kc.l2, labels), "\n\n")

#### TADPole clustering
kc.tadp <- dtwclust(data, type = "tadpole", k = 4,
                    dc = 1.5, control = ctrl)
cat("Rand index for TADPole:", randIndex(kc.tadp, labels), "\n\n")
plot(kc.tadp)

# Modify plot
plot(kc.tadp, clus = 3:4, labs.arg = list(title = "TADPole, clusters 3 and 4",
                                          x = "time", y = "series"))

#### Registering a custom distance with the 'proxy' package and using it
# Normalized DTW distance
ndtw <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = symmetric2,
              distance.only = TRUE, ...)$normalizedDistance
}

# Registering the function with 'proxy'
if (!pr_DB$entry_exists("nDTW"))
     proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "Normalized DTW with L1 norm")

# Subset of (original) data for speed
kc.ndtw <- dtwclust(CharTraj[31:40], distance = "nDTW",
                    control = ctrl, seed = 8319)
cat("Rand index for nDTW (subset):",
    randIndex(kc.ndtw, CharTrajLabels[31:40]), "\n\n")
plot(kc.ndtw)

#### Hierarchical clustering based on shabe-based distance (different lengths)
hc.sbd <- dtwclust(CharTraj, type = "hierarchical",
                   method = c("average", "single"),
                   k = 20, distance = "sbd", control = ctrl)
plot(hc.sbd[[1]])

# Update parameters and re-use distmat
hc.sbd2 <- update(hc.sbd[[1]], method = "complete", distmat = hc.sbd[[1]]@distmat)
plot(hc.sbd2, type = "series")

\dontrun{
     #### Saving and modifying the ggplot object with custom time
     t <- seq(Sys.Date(), len = 180, by = "day")
     gkc <- plot(kc.l2, time = t, plot = FALSE)

     require(scales)
     gkc + scale_x_date(labels = date_format("%b-%Y"),
                        breaks = date_breaks("2 months"))

     #### Using parallel computation to optimize several random repetitions
     #### and distance matrix calculation
     require(doParallel)

     # Create parallel workers
     cl <- makeCluster(detectCores())
     invisible(clusterEvalQ(cl, library(dtwclust)))
     registerDoParallel(cl)

     ctrl <- new("dtwclustControl")
     ctrl@trace <- TRUE

     ## Use full DTW and PAM
     kc.dtw <- dtwclust(CharTraj, k = 20, seed = 3251, control = ctrl)

     ## Use full DTW with DBA centroids
     kc.dba <- dtwclust(CharTraj, k = 20, centroid = "dba", seed = 3251, control = ctrl)

     ## Use constrained DTW with original series of different lengths
     ctrl@window.size <- 20L
     kc.cdtw <- dtwclust(CharTraj, k = 20, seed = 3251, control = ctrl)

     ## This uses the "nDTW" function registered in another example above
     # For reference, this took around 2.25 minutes with 8 cores (all 8 repetitions).
     kc.ndtw.list <- dtwclust(CharTraj, k = 20, distance = "nDTW", centroid = "dba",
                              preproc = zscore, seed = 8319,
                              control = list(window.size = 10L, nrep = 8L))

     # Stop parallel workers
     stopCluster(cl)

     # Return to sequential computations
     registerDoSEQ()

     # See Rand Index for each repetition
     sapply(kc.ndtw.list, randIndex, y = CharTrajLabels)
}
