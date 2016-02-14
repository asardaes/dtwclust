<!-- README.md is generated from README.Rmd. Please edit that file -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dtwclust)](http://cran.r-project.org/web/packages/dtwclust) [![Downloads](http://cranlogs.r-pkg.org/badges/dtwclust)](http://cran.rstudio.com/package=dtwclust)

Time Series Clustering With Dynamic Time Warping Distance (DTW)
===============================================================

This package attempts to consolidate some of the recent techniques related to time series clustering under DTW and implement them in `R`. Most of these algorithms make use of traditional clustering techniques (partitional and hierarchical clustering) but change the distance definition. In this case, the distance between time series is measured with DTW.

DTW is, however, computationally expensive, so several optimization techniques exist. They mostly deal with bounding the DTW distance. These bounds are only defined for time series of equal lengths. Nevertheless, if the length of the time series of interest vary only slightly, reinterpolating them to a common length is probably appropriate.

Additionally, a recently proposed algorithm called k-Shape could serve as an alternative. k-Shape clustering relies on custom distance and centroid definitions, which are unrelated to DTW. The shape extraction algorithm proposed therein is particularly interesting if time series can be z-normalized.

Many of the algorithms and optimizations require that all series have the same length. The ones that don't are usually slow but can still be used.

Please see the included references for more information.

Implementations
---------------

-   Keogh's and Lemire's lower bounds
-   DTW Barycenter Averaging
-   k-Shape clustering
-   TADPole clustering
-   Fuzzy c-means

Examples
--------

``` r
## Load data
data(uciCT)

## Reinterpolate data to equal lengths
datalist <- zscore(CharTraj)
data <- lapply(CharTraj, reinterpolate, newLength = 180)

## Common controls
ctrl <- list(window.size = 20L, trace = TRUE)

## =============================================================================================
## Using DTW with help of lower bounds and PAM centroids
## =============================================================================================

ctrl$pam.precompute <- FALSE

kc.dtwlb <- dtwclust(data = data, k = 20, distance = "dtw_lb",
                     centroid = "pam", seed = 3247, 
                     control = ctrl)
#> Iteration 1: Changes / Distsum = 100 / 1728.288
#> Iteration 2: Changes / Distsum = 18 / 1445.394
#> Iteration 3: Changes / Distsum = 6 / 1360.645
#> Iteration 4: Changes / Distsum = 0 / 1353.769
#> 
#>  Elapsed time is 9.401 seconds.

plot(kc.dtwlb)
```

![](README-examples-1.png)<!-- -->

``` r

ctrl$pam.precompute <- TRUE

## =============================================================================================
## Hierarchical clustering based on shape-based distance
## =============================================================================================

hc.sbd <- dtwclust(datalist, type = "hierarchical",
                   k = 19:21, distance = "sbd",
                   method = "all",
                   control = ctrl)
#> 
#>  Calculating distance matrix...
#> 
#>  Performing hierarchical clustering...
#> 
#>  Elapsed time is 0.65 seconds.

cat("Rand index for HC+SBD:\n")
#> Rand index for HC+SBD:
print(ri <- sapply(hc.sbd, randIndex, y = CharTrajLabels))
#>       ARI       ARI       ARI       ARI       ARI       ARI       ARI 
#> 0.7497894 0.6880051 0.3099641 0.5202854 0.4496453 0.5472890 0.3737120 
#>       ARI       ARI       ARI       ARI       ARI       ARI       ARI 
#> 0.4252906 0.7405012 0.7034438 0.3428238 0.5302644 0.5125830 0.5746408 
#>       ARI       ARI       ARI       ARI       ARI       ARI       ARI 
#> 0.4182468 0.4196458 0.7442112 0.7050215 0.3501627 0.5223116 0.5073598 
#>       ARI       ARI       ARI 
#> 0.5698466 0.4517139 0.4462202

plot(hc.sbd[[which.max(ri)]])
```

![](README-examples-2.png)<!-- -->

``` r

## =============================================================================================
## TADPole clustering
## =============================================================================================

kc.tadp <- dtwclust(data, type = "tadpole", k = 20,
                    dc = 1.5, control = ctrl)
#> 
#> Entering TADPole...
#> 
#> TADPole completed, pruning percentage = 86.7%
#> 
#>  Elapsed time is 4.588 seconds.

plot(kc.tadp, clus = 1:4)
```

![](README-examples-3.png)<!-- -->

``` r

## =============================================================================================
## Parallel support
## =============================================================================================

require(doParallel)
#> Loading required package: doParallel
#> Loading required package: foreach
#> Loading required package: iterators
cl <- makeCluster(detectCores(), "FORK")
invisible(clusterEvalQ(cl, library(dtwclust)))
registerDoParallel(cl)

## Registering a custom distance with proxy and using it (normalized DTW)
ndtw <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = symmetric2,
              distance.only = TRUE, ...)$normalizedDistance
}

## Registering the function with 'proxy'
proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Normalized DTW with L1 norm")

## Data with different lengths
kc.ndtw <- dtwclust(datalist, k = 20,
                    distance = "nDTW", centroid = "pam",
                    seed = 159, control = new("dtwclustControl", nrep = 8L))

sapply(kc.ndtw, randIndex, y = CharTrajLabels)
#>       ARI       ARI       ARI       ARI       ARI       ARI       ARI 
#> 0.6441685 0.5214422 0.4195909 0.4817210 0.5203554 0.5008337 0.5191187 
#>       ARI 
#> 0.5285412

## DBA centroids
kc <- dtwclust(datalist, k = 20,
               distance = "nDTW", centroid = "dba",
               seed = 9421, control = list(trace = TRUE))
#> Series have different lengths. Please confirm that the provided distance function supports this.
#> Iteration 1: Changes / Distsum = 100 / 5.057696
#> Iteration 2: Changes / Distsum = 2 / 3.594286
#> Iteration 3: Changes / Distsum = 1 / 3.550964
#> Iteration 4: Changes / Distsum = 0 / 3.531171
#> 
#>  Elapsed time is 20.315 seconds.

# Modifying some plot parameters
plot(kc, labs.arg = list(title = "DBA Centroids", x = "time", y = "series"))
```

![](README-examples-4.png)<!-- -->

``` r

stopCluster(cl)
registerDoSEQ()

## =============================================================================================
## Fuzzy clustering (autocorrelation-based)
## =============================================================================================

# Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat) {
     lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# Fuzzy distance to be used (squared Euclidean)
fdist <- function(x, y, ...) { sum((x - y)^2) }

# Register it with proxy
if (!pr_DB$entry_exists("SquaredL2"))
     pr_DB$set_entry(FUN = fdist, names = "SquaredL2",
                     loop = TRUE, type = "metric", distance = TRUE)

# Fuzzy c-means
fc <- dtwclust(datalist[1:25], type = "fuzzy", k = 5,
               preproc = acf_fun, distance = "SquaredL2",
               seed = 123)

fc
#> dtwclust(data = datalist[1:25], type = "fuzzy", k = 5, distance = "SquaredL2", 
#>     preproc = acf_fun, seed = 123)
#> 
#> fuzzy clustering with 5 clusters
#> Using SquaredL2 distance
#> Using acf_fun preprocessing
#> 
#> Time required for analysis:
#>    user  system elapsed 
#>   0.235   0.000   0.235 
#> 
#> Head of fuzzy memberships:
#> 
#>       cluster_1   cluster_2  cluster_3  cluster_4 cluster_5
#> A.V1 0.04466135 0.015281712 0.06143803 0.02852775 0.8500912
#> A.V2 0.02639208 0.007424337 0.03720070 0.01505663 0.9139263
#> A.V3 0.03875895 0.007279209 0.03754466 0.01380372 0.9026135
#> A.V4 0.09141695 0.193520982 0.10580764 0.19876731 0.4104871
#> A.V5 0.09256185 0.162738110 0.11853618 0.17478274 0.4513811
#> B.V1 0.40060488 0.034563465 0.34977154 0.07884391 0.1362162
```

Dependencies
------------

-   Partitional procedures are inspired by the `flexclust` package.
-   Hierarchical procedures use the native `hclust` function.
-   Cross-distances make use of the `proxy` package.
-   The core DTW calculations are done by the `dtw` package.
-   Plotting is done with the `ggplot2` package.
-   Parallel computation depends on the `foreach` package.
-   Random streams for repetitions of partitional procedures use the `doRNG` package.
