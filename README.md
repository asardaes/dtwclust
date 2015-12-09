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
data <- zscore(data)

#### Using DTW with help of lower bounds and PAM centroids
kc <- dtwclust(data = data, k = 20, distance = "dtw_lb",
               window.size = 20, centroid = "pam",
               seed = 3247, trace = TRUE)
#>      1 Changes / Distsum : 100 / 1179.209 
#>      2 Changes / Distsum : 13 / 857.7509 
#>      3 Changes / Distsum : 2 / 857.9886 
#>      4 Changes / Distsum : 0 / 857.9886 
#> 
#>  Elapsed time is 5.908 seconds.

plot(kc)
```

![](README-examples-1.png)

``` r

#### Registering a custom distance with proxy and using it (normalized DTW)
ndtw <- function(x, y, ...) {
  dtw::dtw(x, y, step.pattern = symmetric2,
           distance.only = TRUE, ...)$normalizedDistance
}

# Registering the function with 'proxy'
proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Normalized DTW with L1 norm")

# Subset of data for speed (but using different lengths)
kc.ndtw <- dtwclust(datalist[21:40], k = 4, distance = "nDTW",
                    trace = TRUE, seed = 8319)
#>      1 Changes / Distsum : 20 / 2.520941 
#>      2 Changes / Distsum : 6 / 0.8101177 
#>      3 Changes / Distsum : 0 / 0.7401798 
#> 
#>  Elapsed time is 2.228 seconds.

# Modifying some plot parameters
plot(kc.ndtw, labs.arg = list(title = "nDTW clustering", x = "time", y = "series"))
```

![](README-examples-2.png)

``` r

#### Hierarchical clustering based on shape-based distance
hc.sbd <- dtwclust(datalist, type = "hierarchical",
                   distance = "sbd", trace = TRUE)
#> 
#>  Calculating distance matrix...
#> 
#>  Performing hierarchical clustering...
#> 
#>  Elapsed time is 0.814 seconds.
cl.sbd <- cutree(hc.sbd, 20)
cat("Rand index for HC+SBD:", randIndex(cl.sbd, CharTrajLabels), "\n\n")
#> Rand index for HC+SBD: 0.512583
plot(hc.sbd)
```

![](README-examples-3.png)

``` r

#### TADPole clustering
kc.tadp <- dtwclust(data, type = "tadpole", k = 20,
                    window.size = 20, dc = 1.5, 
                    trace = TRUE)
#> 
#> Entering TADPole...
#> 
#> TADPole completed, pruning percentage = 78.5%
#> 
#>  Elapsed time is 8.418 seconds.

plot(kc.tadp, clus = 1:4)
```

![](README-examples-4.png)

Dependencies
------------

-   Partitional procedures are implemented by leveraging the `flexclust` package.
-   Hierarchical procedures use the native `hclust` function.
-   Cross-distances make use of the `proxy` package.
-   The core DTW calculations are done by the `dtw` package.
-   Plotting is done with the `ggplot2` package.
-   Parallel computation depends on the `foreach` package.
