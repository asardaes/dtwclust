<!-- README.md is generated from README.Rmd. Please edit that file -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dtwclust)](https://cran.r-project.org/package=dtwclust) [![Downloads](http://cranlogs.r-pkg.org/badges/dtwclust)](https://cran.r-project.org/package=dtwclust)

Time Series Clustering Along with Optimizations for the Dynamic Time Warping (DTW) Distance
===========================================================================================

Time series clustering with a wide variety of strategies and a series of optimizations specific to the Dynamic Time Warping (DTW) distance and its corresponding lower bounds (LBs). There are implementations of both traditional clustering algorithms, and more recent procedures such as k-Shape and TADPole clustering. Functionality can be easily extended with custom distance measures and centroid definitions.

Many of the algorithms implemented in this package are specifically tailored to time series and DTW, hence its name. However, the main clustering function is flexible so that one can test many different clustering approaches, using either the time series directly, or by applying suitable transformations and then clustering in the resulting space.

DTW is a dynamic programming algorithm that tries to find the optimum warping path between two series. Over the years, several variations have appeared in order to make the procedure faster or more efficient. Please refer to the included references for more information, especially Giorgino (2009), which is a good practical introduction.

Most optimizations require equal dimensionality, which means time series should have equal length. DTW itself does not require this, but it is relatively expensive to compute. Other distance definitions may be used, or series could be reinterpolated to a matching length (Ratanamahatana and Keogh, 2004).

Implementations
---------------

-   Partitional, hierarchical and fuzzy clustering
-   k-Shape clustering
-   TADPole clustering
-   DTW Barycenter Averaging
-   Keogh's and Lemire's DTW lower bounds

Installation
------------

The latest version from CRAN can be installed with `install.packages("dtwclust")`.

If you want to test the latest version from github, first install the [prerequisites for R package development](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) and then type `devtools::install_github("asardaes/dtwclust")`. If you want the vignette to be installed, set the `build_vignettes` parameter to `TRUE`, it will take a couple of minutes.

If you're wondering about which version to install, take a look at the [NEWS](NEWS.md) file, it contains the changelog and I try to keep it updated.

Examples
--------

``` r
## Load data
data(uciCT)

## Reinterpolate data to equal length
datalist <- zscore(CharTraj)
data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))

## Common controls
ctrl <- new("dtwclustControl", window.size = 20L, trace = TRUE)
```

### Partitional

``` r
## =============================================================================================
## Partitional clustering using DTW with help of lower bounds and PAM centroids
## =============================================================================================

ctrl@pam.precompute <- FALSE

kc.dtwlb <- dtwclust(data = data, k = 20, distance = "dtw_lb",
                     centroid = "pam", seed = 3247, 
                     control = ctrl)
#> Iteration 1: Changes / Distsum = 100 / 2055.902
#> Iteration 2: Changes / Distsum = 11 / 1743.629
#> Iteration 3: Changes / Distsum = 5 / 1702.148
#> Iteration 4: Changes / Distsum = 3 / 1702.929
#> Iteration 5: Changes / Distsum = 2 / 1690.314
#> Iteration 6: Changes / Distsum = 0 / 1690.314
#> 
#>  Elapsed time is 5.698 seconds.

plot(kc.dtwlb)
```

![](README-partitional-1.png)

### Hierarchical

``` r
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
#>  Elapsed time is 0.807 seconds.

cat("CVIs for HC+SBD:\n")
#> CVIs for HC+SBD:
print(cvis <- sapply(hc.sbd, cvi, b = CharTrajLabels))
#>               [,1]        [,2]       [,3]        [,4]        [,5]
#> ARI     0.74978939  0.68800515  0.3099641  0.52028537  0.44964526
#> RI      0.97818182  0.97030303  0.8824242  0.94525253  0.92686869
#> J       0.61428571  0.54205607  0.2145749  0.37557604  0.31698113
#> FM      0.76614889  0.71634666  0.4249454  0.57846535  0.53232837
#> VI      0.19822365  0.22966355  0.5661339  0.35605205  0.39562319
#> Sil     0.57988078  0.60407120  0.3888051  0.55712005  0.54459348
#> SF      0.36643161  0.38281683  0.4701455  0.41542295  0.44875986
#> CH     27.05423386 26.54335013 16.9670382 24.03940574 22.49842911
#> DB      0.63583376  0.64043681  0.5217145  0.61107493  0.49464680
#> DBstar  1.61936383  1.11335876  0.9097630  0.78382455  0.69307237
#> D       0.13801520  0.21341626  0.1308229  0.18439480  0.17588348
#> COP     0.06956278  0.07242208  0.1244863  0.08522869  0.09168347
#>               [,6]       [,7]        [,8]        [,9]       [,10]
#> ARI     0.54728896  0.3737120  0.42529063  0.74050121  0.70344375
#> RI      0.94909091  0.9074747  0.92181818  0.97838384  0.97232323
#> J       0.40000000  0.2588997  0.29764065  0.60223048  0.55948553
#> FM      0.60308485  0.4705882  0.51100510  0.75369221  0.72880580
#> VI      0.34983127  0.4936223  0.43496098  0.22265488  0.21511035
#> Sil     0.55292775  0.4498256  0.50365682  0.58811104  0.61748556
#> SF      0.42993257  0.4555955  0.45353987  0.36348668  0.38500224
#> CH     22.43963739 18.1705294 19.68318121 27.29320621 25.53616172
#> DB      0.62649917  0.5316269  0.46641999  0.65336293  0.61402700
#> DBstar  0.77189375  0.7499572  0.65255819  1.68466461  1.09582372
#> D       0.16176333  0.1747119  0.19373597  0.13801520  0.22382404
#> COP     0.08582413  0.1085514  0.09651714  0.06595574  0.06893191
#>             [,11]       [,12]       [,13]       [,14]      [,15]
#> ARI     0.3428238  0.53026441  0.51258299  0.57464085  0.4182468
#> RI      0.8965657  0.94848485  0.94202020  0.95414141  0.9216162
#> J       0.2369598  0.38405797  0.36923077  0.42531646  0.2919708
#> FM      0.4479318  0.58214036  0.57759590  0.62350648  0.5019646
#> VI      0.5185772  0.34855573  0.34677938  0.31972827  0.4460655
#> Sil     0.4094015  0.53663563  0.56617496  0.57151588  0.4741491
#> SF      0.4659191  0.40333496  0.44479800  0.43608567  0.4518092
#> CH     17.1076283 23.50969213 23.28285129 23.03569943 19.6161099
#> DB      0.5085174  0.72846310  0.46973394  0.55074053  0.5464752
#> DBstar  0.9083771  0.95600061  0.61732267  0.64117537  0.7977064
#> D       0.1308103  0.17967703  0.23163554  0.16176333  0.1747119
#> COP     0.1154545  0.08337437  0.08273567  0.07841355  0.1000334
#>              [,16]       [,17]       [,18]      [,19]       [,20]
#> ARI     0.41964584  0.74421123  0.70502147  0.3501627  0.52231158
#> RI      0.92202020  0.97878788  0.97252525  0.9030303  0.94787879
#> J       0.29304029  0.60674157  0.56129032  0.2417062  0.37681159
#> FM      0.50295569  0.75697629  0.73008778  0.4469178  0.57346741
#> VI      0.44257501  0.21436184  0.20908975  0.5179086  0.35832449
#> Sil     0.50197422  0.59978230  0.62335691  0.4332002  0.53741732
#> SF      0.45568653  0.37749889  0.41900332  0.4622818  0.40959973
#> CH     18.93844087 27.17714225 25.39711110 15.8362781 22.56431559
#> DB      0.45639546  0.61387828  0.55601385  0.4538754  0.69863982
#> DBstar  0.66318416  1.22942083  0.70446747  0.8258729  0.93239293
#> D       0.19373597  0.15250536  0.24447014  0.1516362  0.18010010
#> COP     0.09445372  0.06193398  0.06600764  0.1102822  0.08171361
#>              [,21]       [,22]       [,23]       [,24]
#> ARI     0.50735985  0.56984658  0.45171395  0.44622023
#> RI      0.94222222  0.95434343  0.93070707  0.92929293
#> J       0.36444444  0.42051282  0.31809145  0.31372549
#> FM      0.56993940  0.61634974  0.52579262  0.52186246
#> VI      0.35439340  0.32734230  0.40643794  0.40772660
#> Sil     0.56463229  0.56980939  0.51550447  0.52593013
#> SF      0.44692615  0.43820045  0.45165102  0.45082644
#> CH     22.54062825 22.33204629 19.53887943 19.38570054
#> DB      0.46193786  0.53595538  0.49555667  0.46748541
#> DBstar  0.62239300  0.62757437  0.78276218  0.68310239
#> D       0.23163554  0.16176333  0.17471193  0.24428149
#> COP     0.08067224  0.07635012  0.09213935  0.08780314

plot(hc.sbd[[which.min(cvis["VI", ])]])
```

![](README-hierarchical-1.png)

### TADPole

``` r
## =============================================================================================
## TADPole clustering
## =============================================================================================

kc.tadp <- dtwclust(data, type = "tadpole", k = 20,
                    dc = 1.5, control = ctrl)
#> 
#> Entering TADPole...
#> 
#> TADPole completed, pruning percentage = 86.5%
#> 
#>  Elapsed time is 1.584 seconds.

plot(kc.tadp, clus = 1:4)
```

![](README-tadpole-1.png)

### Parallel support

``` r
## =============================================================================================
## Parallel support
## =============================================================================================

require(doParallel)
#> Loading required package: doParallel
#> Loading required package: foreach
#> Loading required package: iterators
cl <- makeCluster(detectCores(), "FORK")
registerDoParallel(cl)

## Creating a custom distance (normalized DTW)
ndtw <- function(x, y, ...) {
     dtw::dtw(x, y, step.pattern = symmetric2,
              distance.only = TRUE, ...)$normalizedDistance
}

## Registering the function with 'proxy'
proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Normalized DTW with L1 norm")

## Data with different length
kc.ndtw <- dtwclust(datalist, k = 20,
                    distance = "nDTW", centroid = "pam",
                    seed = 159, control = new("dtwclustControl", nrep = 8L))

sapply(kc.ndtw, cvi, b = CharTrajLabels, type = "VI")
#>        VI        VI        VI        VI        VI        VI        VI 
#> 0.2755874 0.4448495 0.5736902 0.5047190 0.4465854 0.3929220 0.3808492 
#>        VI 
#> 0.4439003

## DBA centroids
kc <- dtwclust(datalist, k = 20,
               distance = "dtw_basic", centroid = "dba",
               seed = 9421, control = list(trace = TRUE, window.size = 20L),
               normalize = TRUE)
#> Iteration 1: Changes / Distsum = 100 / 6.749242
#> Iteration 2: Changes / Distsum = 2 / 4.229023
#> Iteration 3: Changes / Distsum = 0 / 4.180198
#> 
#>  Elapsed time is 5.417 seconds.

## Modifying some plot parameters
plot(kc, labs.arg = list(title = "DBA Centroids", x = "time", y = "series"))
```

![](README-parallel-1.png)

``` r

## Returning to sequential calculations
stopCluster(cl)
registerDoSEQ()
```

### Fuzzy

``` r
## =============================================================================================
## Fuzzy clustering (autocorrelation-based)
## =============================================================================================

# Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat, ...) {
     lapply(dat, function(x) as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf))
}

# Fuzzy c-means
fc <- dtwclust(datalist[1:25], type = "fuzzy", k = 5,
               preproc = acf_fun, distance = "L2",
               seed = 123)

fc
#> dtwclust(data = datalist[1:25], type = "fuzzy", k = 5, distance = "L2", 
#>     preproc = acf_fun, seed = 123)
#> 
#> fuzzy clustering with 5 clusters
#> Using L2 distance
#> Using acf_fun preprocessing
#> 
#> Time required for analysis:
#>    user  system elapsed 
#>   0.156   0.000   0.157 
#> 
#> Head of fuzzy memberships:
#> 
#>       cluster_1   cluster_2  cluster_3  cluster_4 cluster_5
#> A.V1 0.04517433 0.015248385 0.06048626 0.02847461 0.8506164
#> A.V2 0.02648172 0.007341668 0.03623648 0.01489308 0.9150471
#> A.V3 0.03920172 0.007216578 0.03668630 0.01368817 0.9032072
#> A.V4 0.09258928 0.193779128 0.10495491 0.19932425 0.4093524
#> A.V5 0.09366124 0.162965470 0.11758524 0.17523731 0.4505507
#> B.V1 0.39400450 0.034717343 0.35507763 0.07914583 0.1370547
```

Dependencies
------------

-   Partitional procedures are inspired by the `flexclust` package.
-   Hierarchical procedures use the native `hclust` function.
-   Cross-distances make use of the `proxy` package.
-   The core DTW calculations are done by the `dtw` package.
-   Plotting is done with the `ggplot2` package.
-   Parallel computation depends on the `foreach` package.
