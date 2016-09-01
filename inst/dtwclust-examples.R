## NOTE: All examples should work, some are simply set to not run automatically

# Load data
data(uciCT)

# Controls
ctrl <- new("dtwclustControl", trace = TRUE)

# ====================================================================================
# Simple partitional clustering with Euclidean distance and PAM centroids
# ====================================================================================

# Reinterpolate to same length and coerce as matrix
data <- t(sapply(CharTraj, reinterpolate, newLength = 180))

# Subset for speed
data <- data[1:20, ]
labels <- CharTrajLabels[1:20]

# Testing several values of k
kc.l2 <- dtwclust(data, k = 2:4,
                  distance = "L2", centroid = "pam",
                  seed = 3247, control = ctrl)

cat("Rand indices for L2+PAM:\n")
sapply(kc.l2, randIndex, y = labels)

plot(kc.l2[[3L]])

# ====================================================================================
# Hierarchical clustering with Euclidean distance
# ====================================================================================

# Re-use the distance matrix from the previous example
hc.l2 <- update(kc.l2[[3L]], k = 4L,
                type = "hierarchical", method = "all",
                distmat = kc.l2[[3L]]@distmat)

cat("VI for L2+HC:\n")
sapply(hc.l2, cvi, b = labels, type = "VI")

# ====================================================================================
# Registering a custom distance with the 'proxy' package and using it
# ====================================================================================

# Normalized DTW distance
ndtw <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = symmetric2,
              distance.only = TRUE, ...)$normalizedDistance
}

# Registering the function with 'proxy'
if (!pr_DB$entry_exists("nDTW"))
     proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                            loop = TRUE, type = "metric", distance = TRUE,
                            description = "Normalized DTW")

# Subset of (original) data for speed
kc.ndtw <- dtwclust(CharTraj[31:39], k = 2L,
                    distance = "nDTW",
                    control = ctrl, seed = 8319)

# Where would series 40 go?
predict(kc.ndtw, newdata = CharTraj[40])

# ====================================================================================
# Hierarchical clustering based on shabe-based distance
# ====================================================================================

# Original data
hc.sbd <- dtwclust(CharTraj, type = "hierarchical",
                   k = 20, method = "all",
                   distance = "sbd", control = ctrl)

# Plot dendrogram
plot(hc.sbd[[ which.max(sapply(hc.sbd, randIndex, y = CharTrajLabels)) ]])

# ====================================================================================
# Autocorrelation-based fuzzy clustering (see D'Urso and Maharaj 2009)
# ====================================================================================

# Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat) {
     lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# Fuzzy c-means
fc <- dtwclust(CharTraj[1:25], type = "fuzzy", k = 5,
               preproc = acf_fun, distance = "L2",
               seed = 123, control = ctrl)

# Fuzzy memberships
print(fc@fcluster)

# Plot crisp partition in the original space
plot(fc, data = CharTraj[1:25], type = "series")

# Membership for new series
# (It's important that 'preproc' was specified in the original call)
predict(fc, newdata = CharTraj[26:28])

\dontrun{
# ====================================================================================
# TADPole clustering
# ====================================================================================

ctrl@window.size <- 20L

kc.tadp <- dtwclust(data, type = "tadpole", k = 4L,
                    dc = 1.5, control = ctrl)

cat("Rand index for TADPole:", randIndex(kc.tadp, labels), "\n\n")

# Modify plot
plot(kc.tadp, clus = 3:4,
     labs.arg = list(title = "TADPole, clusters 3 and 4",
                     x = "time", y = "series"))

# ====================================================================================
# Saving and modifying the ggplot object with custom time labels
# ====================================================================================

t <- seq(Sys.Date(), len = 180, by = "day")
gkc <- plot(kc.tadp, time = t, plot = FALSE)

require(scales)
gkc + scale_x_date(labels = date_format("%b-%Y"),
                   breaks = date_breaks("2 months"))

# ====================================================================================
# Specifying a centroid function for prototype extraction in hierarchical clustering
# ====================================================================================

# Seed is due to possible randomness in shape_extraction when selecting a basis series
hc <- dtwclust(CharTraj, type = "hierarchical",
               k = 20, distance = "sbd",
               preproc = zscore, centroid = shape_extraction,
               seed = 320)

plot(hc, type = "sc")

# ====================================================================================
# Multivariate time series
# ====================================================================================

# Dummy multivariate series
mv <- lapply(seq(1L, 100L, 5L), function(x) cbind(CharTraj[[x]], CharTraj[[x+1L]]))

# "dist.method" is for the dtw function, and it affects the values of the local
# cost matrix in the case of multivariate series
mvc <- dtwclust(mv, k = 4L, dist.method = "L1", seed = 390)

# Note how the "dimensions" of each series are appended one after the other in the plot
plot(mvc)

# ====================================================================================
# Using parallel computation to optimize several random repetitions
# and distance matrix calculation
# ====================================================================================
require(doParallel)

# Create parallel workers
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, library(dtwclust)))
registerDoParallel(cl)

## Use constrained DTW and PAM (~20 seconds with 8 cores)
kc.dtw <- dtwclust(CharTraj, k = 20, seed = 3251, control = ctrl)

## Use constrained DTW with DBA centroids (~20 seconds with 8 cores)
kc.dba <- dtwclust(CharTraj, k = 20, centroid = "dba", seed = 3251, control = ctrl)

## This uses the "nDTW" function registered in another example above
# For reference, this took around 2.25 minutes with 8 cores (all 8 repetitions).
kc.ndtw.reps <- dtwclust(CharTraj, k = 20, distance = "nDTW", centroid = "dba",
                         preproc = zscore, seed = 8319,
                         control = list(window.size = 10L, nrep = 8L))

# See VI index for each repetition
sapply(kc.ndtw.reps, cvi, b = CharTrajLabels, type = "VI")

# Stop parallel workers
stopCluster(cl)

# Return to sequential computations
registerDoSEQ()
}
