<!-- README.md is generated from README.Rmd. Please edit that file -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dtwclust)](https://cran.r-project.org/package=dtwclust) [![Downloads](http://cranlogs.r-pkg.org/badges/dtwclust)](https://cran.r-project.org/package=dtwclust)

Time Series Clustering Along with Optimizations for the Dynamic Time Warping Distance (DTW)
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

Examples
--------

``` r
## Load data
data(uciCT)

## Reinterpolate data to equal length
datalist <- zscore(CharTraj)
data <- lapply(CharTraj, reinterpolate, newLength = 180)

## Common controls
ctrl <- new("dtwclustControl", window.size = 20L, trace = TRUE)
```

``` r
## =============================================================================================
## Partitional clustering using DTW with help of lower bounds and PAM centroids
## =============================================================================================

ctrl@pam.precompute <- FALSE

kc.dtwlb <- dtwclust(data = data, k = 20, distance = "dtw_lb",
                     centroid = "pam", seed = 3247, 
                     control = ctrl)
#> Iteration 1: Changes / Distsum = 100 / 1747.417
#> Iteration 2: Changes / Distsum = 18 / 1417.733
#> Iteration 3: Changes / Distsum = 13 / 1349.521
#> Iteration 4: Changes / Distsum = 2 / 1311.201
#> Iteration 5: Changes / Distsum = 0 / 1311.201
#> 
#>  Elapsed time is 9.528 seconds.

plot(kc.dtwlb)
```

![](README-partitional-1.png)

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
#>  Elapsed time is 0.649 seconds.

cat("CVIs for HC+SBD:\n")
#> CVIs for HC+SBD:
print(cvis <- sapply(hc.sbd, cvi, b = CharTrajLabels))
#>              [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> ARI     0.7497894  0.6880051  0.3099641  0.5202854  0.4496453  0.5472890
#> RI      0.9781818  0.9703030  0.8824242  0.9452525  0.9268687  0.9490909
#> J       0.6142857  0.5420561  0.2145749  0.3755760  0.3169811  0.4000000
#> FM      0.7661489  0.7163467  0.4249454  0.5784654  0.5323284  0.6030849
#> VI      0.1982237  0.2296635  0.5661339  0.3560521  0.3956232  0.3498313
#> Sil     0.5798808  0.6040712  0.3888051  0.5571201  0.5445935  0.5529278
#> SF      0.3664316  0.3828168  0.4701455  0.4154229  0.4487599  0.4299326
#> CH     27.0542339 26.5433501 16.9670382 24.0394057 22.4984291 22.4396374
#> DB      0.6358338  0.6404368  0.5217145  0.6110749  0.4946468  0.6264992
#> DBstar  1.6193638  1.1133588  0.9097630  0.7838245  0.6930724  0.7718938
#> D       0.1380152  0.2134163  0.1308229  0.1843948  0.1758835  0.1617633
#>              [,7]       [,8]       [,9]      [,10]      [,11]      [,12]
#> ARI     0.3737120  0.4252906  0.7405012  0.7034438  0.3428238  0.5302644
#> RI      0.9074747  0.9218182  0.9783838  0.9723232  0.8965657  0.9484848
#> J       0.2588997  0.2976407  0.6022305  0.5594855  0.2369598  0.3840580
#> FM      0.4705882  0.5110051  0.7536922  0.7288058  0.4479318  0.5821404
#> VI      0.4936223  0.4349610  0.2226549  0.2151104  0.5185772  0.3485557
#> Sil     0.4498256  0.5036568  0.5881110  0.6174856  0.4094015  0.5366356
#> SF      0.4555955  0.4535399  0.3634867  0.3850022  0.4659191  0.4033350
#> CH     18.1705294 19.6831812 27.2932062 25.5361617 17.1076283 23.5096921
#> DB      0.5316269  0.4664200  0.6533629  0.6140270  0.5085174  0.7284631
#> DBstar  0.7499572  0.6525582  1.6846646  1.0958237  0.9083771  0.9560006
#> D       0.1747119  0.1937360  0.1380152  0.2238240  0.1308103  0.1796770
#>             [,13]      [,14]      [,15]      [,16]      [,17]      [,18]
#> ARI     0.5125830  0.5746408  0.4182468  0.4196458  0.7442112  0.7050215
#> RI      0.9420202  0.9541414  0.9216162  0.9220202  0.9787879  0.9725253
#> J       0.3692308  0.4253165  0.2919708  0.2930403  0.6067416  0.5612903
#> FM      0.5775959  0.6235065  0.5019646  0.5029557  0.7569763  0.7300878
#> VI      0.3467794  0.3197283  0.4460655  0.4425750  0.2143618  0.2090898
#> Sil     0.5661750  0.5715159  0.4741491  0.5019742  0.5997823  0.6233569
#> SF      0.4447980  0.4360857  0.4518092  0.4556865  0.3774989  0.4190033
#> CH     23.2828513 23.0356994 19.6161099 18.9384409 27.1771422 25.3971111
#> DB      0.4697339  0.5507405  0.5464752  0.4563955  0.6138783  0.5560139
#> DBstar  0.6173227  0.6411754  0.7977064  0.6631842  1.2294208  0.7044675
#> D       0.2316355  0.1617633  0.1747119  0.1937360  0.1525054  0.2444701
#>             [,19]      [,20]      [,21]      [,22]      [,23]      [,24]
#> ARI     0.3501627  0.5223116  0.5073598  0.5698466  0.4517139  0.4462202
#> RI      0.9030303  0.9478788  0.9422222  0.9543434  0.9307071  0.9292929
#> J       0.2417062  0.3768116  0.3644444  0.4205128  0.3180915  0.3137255
#> FM      0.4469178  0.5734674  0.5699394  0.6163497  0.5257926  0.5218625
#> VI      0.5179086  0.3583245  0.3543934  0.3273423  0.4064379  0.4077266
#> Sil     0.4332002  0.5374173  0.5646323  0.5698094  0.5155045  0.5259301
#> SF      0.4622818  0.4095997  0.4469262  0.4382004  0.4516510  0.4508264
#> CH     15.8362781 22.5643156 22.5406282 22.3320463 19.5388794 19.3857005
#> DB      0.4538754  0.6986398  0.4619379  0.5359554  0.4955567  0.4674854
#> DBstar  0.8258729  0.9323929  0.6223930  0.6275744  0.7827622  0.6831024
#> D       0.1516362  0.1801001  0.2316355  0.1617633  0.1747119  0.2442815

plot(hc.sbd[[which.min(cvis["VI", ])]])
```

![](README-hierarchical-1.png)

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
#>  Elapsed time is 3.482 seconds.

plot(kc.tadp, clus = 1:4)
```

![](README-tadpole-1.png)

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
               distance = "nDTW", centroid = "dba",
               seed = 9421, control = list(trace = TRUE))
#> Series have different length. Please confirm that the provided distance function supports this.
#> Iteration 1: Changes / Distsum = 100 / 5.057696
#> Iteration 2: Changes / Distsum = 2 / 3.594286
#> Iteration 3: Changes / Distsum = 1 / 3.550964
#> Iteration 4: Changes / Distsum = 0 / 3.531171
#> 
#>  Elapsed time is 17.91 seconds.

## Modifying some plot parameters
plot(kc, labs.arg = list(title = "DBA Centroids", x = "time", y = "series"))
```

![](README-parallel-1.png)

``` r

## Returning to sequential calculations
stopCluster(cl)
registerDoSEQ()
```

``` r
## =============================================================================================
## Fuzzy clustering (autocorrelation-based)
## =============================================================================================

# Calculate autocorrelation up to 50th lag, considering a list of time series as input
acf_fun <- function(dat) {
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
#>   0.144   0.000   0.143 
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
