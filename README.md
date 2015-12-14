<!-- README.md is generated from README.Rmd. Please edit that file -->
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

Examples
--------

``` r
## Load data
data(uciCT)

## Reinterpolate data to equal lengths
datalist <- zscore(CharTraj)
data <- lapply(CharTraj, reinterpolate, newLength = 180)
# data <- zscore(data)

#### Using DTW with help of lower bounds and PAM centroids
kc.dtwlb <- dtwclust(data = data, k = 20, distance = "dtw_lb",
                     window.size = 20, centroid = "pam",
                     seed = 3247, trace = TRUE)
#> Iteration 1: Changes / Distsum =100 / 932.183
#> Iteration 2: Changes / Distsum =18 / 632.6588
#> Iteration 3: Changes / Distsum =2 / 596.8066
#> Iteration 4: Changes / Distsum =0 / 596.5887
#> 
#> 
#>  Elapsed time is 4.322 seconds.

plot(kc.dtwlb)
```

![](README-examples-1.png)

``` r

#### Hierarchical clustering based on shape-based distance
hc.sbd <- dtwclust(datalist, type = "hierarchical",
                   k = 20, distance = "sbd", trace = TRUE)
#> 
#>  Calculating distance matrix...
#> 
#>  Performing hierarchical clustering...
#> 
#>  Elapsed time is 0.615 seconds.

require(flexclust)
cat("Rand index for HC+SBD:", randIndex(hc.sbd@cluster, CharTrajLabels), "\n\n")
#> Rand index for HC+SBD: 0.512583
plot(hc.sbd)
```

![](README-examples-2.png)

``` r

#### TADPole clustering
kc.tadp <- dtwclust(data, type = "tadpole", k = 20,
                    window.size = 20, dc = 1.5, 
                    trace = TRUE)
#> 
#> Entering TADPole...
#> 
#> TADPole completed, pruning percentage = 86.7%
#> 
#>  Elapsed time is 3.962 seconds.

plot(kc.tadp, clus = 1:4)
```

![](README-examples-3.png)

``` r

#### Parallel support
require(doParallel)
#> Loading required package: doParallel
#> Loading required package: iterators
#> Loading required package: parallel
cl <- makeCluster(detectCores(), "FORK")
registerDoParallel(cl)

## Registering a custom distance with proxy and using it (normalized DTW)
ndtw <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = symmetric2,
              distance.only = TRUE, ...)$normalizedDistance
}

# Registering the function with 'proxy'
proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Normalized DTW with L1 norm")

# Data with different lengths
kc.ndtw <- dtwclust(datalist, k = 20,
                    distance = "nDTW", centroid = "pam",
                    trace = TRUE, seed = 159, reps = 8L)
#> Consider setting save.data to FALSE if performing several repetitions.
#> 
#> Tracing will not be available if parallel computing is used.
#> 
#>  Elapsed time is 5.44 seconds.

sapply(kc.ndtw, function(x) randIndex(x@cluster, CharTrajLabels))
#>       ARI       ARI       ARI       ARI       ARI       ARI       ARI 
#> 0.5739085 0.5798654 0.4112420 0.5351719 0.6000514 0.4897014 0.5101194 
#>       ARI 
#> 0.6715448

# DBA centroids
kc <- dtwclust(datalist, k = 20,
               distance = "nDTW", centroid = "dba",
               trace = TRUE, seed = 9421)
#> Iteration 1: Changes / Distsum =100 / 5.162033
#> Iteration 2: Changes / Distsum =3 / 3.739462
#> Iteration 3: Changes / Distsum =2 / 3.687197
#> Iteration 4: Changes / Distsum =0 / 3.631238
#> 
#> 
#>  Elapsed time is 22.432 seconds.

# Modifying some plot parameters
plot(kc, labs.arg = list(title = "DBA Centroids", x = "time", y = "series"))
```

![](README-examples-4.png)

``` r

stopCluster(cl)
registerDoSEQ()
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
