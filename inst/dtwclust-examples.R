#' NOTE: More examples are available in the vignette. Here are just some miscellaneous
#' examples that might come in handy. They should all work, but some don't run
#' automatically.

# Load data
data(uciCT)

# ====================================================================================
# Simple partitional clustering with Euclidean distance and PAM centroids
# ====================================================================================

# Reinterpolate to same length
series <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))

# Subset for speed
series <- series[1:20]
labels <- CharTrajLabels[1:20]

# Making many repetitions
pc.l2 <- dtwclust(series, k = 4L,
                  distance = "L2", centroid = "pam",
                  seed = 3247, control = list(trace = TRUE,
                                              nrep = 10L))

# Cluster validity indices
sapply(pc.l2, cvi, b = labels)

# ====================================================================================
# Hierarchical clustering with Euclidean distance
# ====================================================================================

# Re-use the distance matrix from the previous example (all matrices are the same)
# Use all available linkage methods for function hclust
hc.l2 <- dtwclust(series, type = "hierarchical",
                  k = 4L, method = "all",
                  control = list(trace = TRUE),
                  distmat = pc.l2[[1L]]@distmat)

# Plot the best dendrogram according to variation of information
plot(hc.l2[[which.min(sapply(hc.l2, cvi, b = labels, type = "VI"))]])

# ====================================================================================
# Multivariate time series
# ====================================================================================

# Multivariate series, provided as a list of matrices
mv <- CharTrajMV[1L:20L]

# Using GAK distance
mvc <- dtwclust(mv, k = 4L, distance = "gak", seed = 390)

# Note how the variables of each series are appended one after the other in the plot
plot(mvc)

\dontrun{
# Common controls
ctrl <- new("dtwclustControl", trace = TRUE, window.size = 18L)

# ====================================================================================
# Registering a custom distance with the 'proxy' package and using it
# ====================================================================================

# Normalized asymmetric DTW distance
ndtw <- function(x, y, ...) {
    dtw::dtw(x, y, step.pattern = asymmetric,
             distance.only = TRUE, ...)$normalizedDistance
}

# Registering the function with 'proxy'
if (!pr_DB$entry_exists("nDTW"))
    proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                           loop = TRUE, type = "metric", distance = TRUE,
                           description = "Normalized asymmetric DTW")

# Subset of (original) data for speed
pc.ndtw <- dtwclust(series[-1L], k = 4L,
                    distance = "nDTW",
                    seed = 8319,
                    control = ctrl)

# Which cluster would the first series belong to?
# Notice that newdata is provided as a list
predict(pc.ndtw, newdata = series[1L])

# ====================================================================================
# Custom hierarchical clustering
# ====================================================================================

require(cluster)

hc.diana <- dtwclust(series, type = "h", k = 4L,
                     distance = "L2", method = diana,
                     control = ctrl)

plot(hc.diana, type = "sc")

# ====================================================================================
# TADPole clustering
# ====================================================================================

pc.tadp <- dtwclust(series, type = "tadpole", k = 4L,
                    dc = 1.5, control = ctrl)

# Modify plot, show only clusters 3 and 4
plot(pc.tadp, clus = 3:4,
     labs.arg = list(title = "TADPole, clusters 3 and 4",
                     x = "time", y = "series"))

# Saving and modifying the ggplot object with custom time labels
require(scales)
t <- seq(Sys.Date(), len = length(series[[1L]]), by = "day")
gpc <- plot(pc.tadp, time = t, plot = FALSE)

gpc + scale_x_date(labels = date_format("%b-%Y"),
                   breaks = date_breaks("2 months"))

# ====================================================================================
# Specifying a centroid function for prototype extraction in hierarchical clustering
# ====================================================================================

# Seed is due to possible randomness in shape_extraction when selecting a basis series
hc.sbd <- dtwclust(CharTraj, type = "hierarchical",
                   k = 20L, distance = "sbd",
                   preproc = zscore, centroid = shape_extraction,
                   seed = 320L)

plot(hc.sbd, type = "sc")

# ====================================================================================
# Using parallel computation to optimize several random repetitions
# and distance matrix calculation
# ====================================================================================
require(doParallel)

# Create parallel workers
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, library(dtwclust)))
registerDoParallel(cl)

## Use constrained DTW and PAM (less than a second with 8 cores)
pc.dtw <- dtwclust(CharTraj, k = 20L, seed = 3251, control = ctrl)

## Use constrained DTW with DBA centroids (~3 seconds with 8 cores)
pc.dba <- dtwclust(CharTraj, k = 20L, centroid = "dba", seed = 3251, control = ctrl)

#' Using distance based on global alignment kernels
#' (~35 seconds with 8 cores for all repetitions)
pc.gak <- dtwclust(CharTraj, k = 20L,
                   distance = "gak",
                   centroid = "dba",
                   seed = 8319,
                   control = list(trace = TRUE,
                                  window.size = 18L,
                                  nrep = 8L))

# Stop parallel workers
stopCluster(cl)

# Return to sequential computations. This MUST be done after stopCluster()
registerDoSEQ()
}
