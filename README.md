<!-- README.md is generated from README.Rmd. Please edit that file -->
Time Series Clustering With Dynamic Time Warping Distance (DTW)
===============================================================

This package attempts to consolidate some of the recent techniques related to time series clustering under DTW and implement them in `R`. Most of these algorithms make use of traditional clustering techniques (partitional and hierarchical clustering) but change the distance definition. In this case, the distance between time series is measured with DTW.

DTW is, however, computationally expensive, so several optimization techniques exist. They mostly deal with bounding the DTW distance. These bounds are only defined for time series of equal lengths. Nevertheless, if the length of the time series of interest vary only slightly, reinterpolating them to a common length is probably appropriate.

Additionally, a recently proposed algorithm called k-Shape could serve as an alternative. k-Shape clustering relies on custom distance and centroid definitions, which are unrelated to DTW. The shape extraction algorithm proposed therein is particularly interesting if time series can be normalized.

Many of the algorithms and optimizations require that all series have the same length. The ones that don't are usually slow but can still be used.

Please see the included references for more information.

Example
-------

``` r
# Load data
data(uciCT)

# Reinterpolate data to equal lengths
data <- lapply(CharTraj, reinterpolate, newLength = 205)

kc <- dtwclust(data = data, k = 20, distance = "dtw_lb",
               window.size = 20, centroid = "pam",
               save.data = TRUE,
               seed = 3247, trace = TRUE)
#>      1 Changes / Distsum : 100 / 1157.388 
#>      2 Changes / Distsum : 17 / 826.3151 
#>      3 Changes / Distsum : 4 / 752.9349 
#>      4 Changes / Distsum : 1 / 752.1322 
#>      5 Changes / Distsum : 0 / 752.1322 
#> 
#>  Elapsed time is 5.53 seconds.

plot(kc)
```

![](README-example-1.png)

Dependencies
------------

-   Partitional procedures are implemented by leveraging the `flexclust` package.
-   Hierarchical procedures use the native `hclust` function.
-   Cross-distances make use of the `proxy` package.
-   The core DTW calculations are done by the `dtw` package.
-   Plotting is done with the `ggplot2` package.

Implementations
---------------

-   Keogh's and Lemire's lower bounds
-   DTW Barycenter Averaging
-   k-Shape clustering
-   TADPole clustering
