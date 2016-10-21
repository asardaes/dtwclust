# NEWS

## Version 3.0.0
* Removed deprecated arguments/slots.
    + If you have older `dtwclust` objects saved, try updating them with `attr(dtwclust_object, "centers") <- NULL` if you run into compatibility problems.
    + `DBA` arguments changed order.
* Fixed detection of some symmetric DTW cases.
* Added (conditional) support for more hierarchical procedures. See the examples.
* Added support for functions in package `clue`.
* Bear in mind that the `DTW/SBD` algorithms (and hence the functions that depend on them) might give different results in 32-bit systems.

## Version 2.3.0
* Correction: DTW can be symmetric for series of both equal/different length, although in general this is not necessarily the case, due to asymmetric step patterns or constrained paths.
* Exported the `dtw2` function.
* Revamped `zscore` and `reinterpolate` functions.
* Data frames are now parsed row-wise, like matrices, to maintain consistency with `proxy`.
* Fixed a bug where multivariate series with different length had spurious data added to them.
* Fixed a bug in multivariate `shape_extraction`. Reminder: multivariate shape extraction might not be a good idea.
* Consistency adjustments: all `center(s)` arguments/slots will be removed and replaced with `centroid(s)`.
* Added the `dtw_basic` function, which should be faster due to its limited functionality. It can
also be used with `dtw_lb` and `DBA`. It will now be used by default.
* The `SBD` and `shape_extraction` functions used the `crossprod` function internally, which returned a 1x1 matrix by default, causing some dimension inconsistencies. In the future, `R` will give an error about the inconsistency, so the function has been changed in `dtwclust`, but it resulted in very small numerical differences, which may be enough to alter some clustering results.

## Version 2.2.1
* Added package vignette
* Implemented several cluster validity indices in the new `cvi` function
* The custom argument `force.pairwise` used in some of the `proxy` distances is not necessary anymore
* No longer enforcing `NULL` dimensions for each series (in case multivariate series are to be used), use with caution
* Fixed an error that prevented fuzzy clustering from working when `k = 2`
* Fixed a possible problem in some of the included proxy distances when a data object had a length of 1

## Version 2.1.2
* Added option to specify a function to extract prototypes in hierarchical and TADPole clustering
* Switched algorithm to calculate envelops to Lemire's, which should be slightly faster
* More plot types
* Documentation fixes

## Version 2.1.1
* Minor bug fixes for fuzzy clustering
     + Fixed the `predict` generic.
     + The final values returned in the `fcluster` slot needed one final update during clustering. It should be correct now, but it will vary slightly with respect to what was previously given in v2.1.0.
* Some more examples

## Version 2.1.0
* Added fuzzy c-means clustering
* Memory optimizations for DBA
* Support for several `k` values for multiple runs
* Update documentation

## Version 2.0.0
* Major refactor.
* Many formal parameters from the `dtwclust` function were dropped and implemented in the formal class `dtwclustControl`. For now, they will still be supported through `...` with a message.
* No longer supporting non-proxy distances, in order to be able to take advantage of the included optimizations.
* Dropped inheritance of `flexclust`'s `kccasimple`
     + Many slots and methods were ported
     + Inheriting from `hclust` class now, and all its associated methods
* `shape_extraction` now accepts series with different lengths!
* More parallel support
     + DBA and Shape centroid calculations 
     + DBA itself (probably unnecessary unless you're averaging a lot of series)
     + Custom `proxy` distances directly (except `DTW2`)
* Several hierarchical procedures can be made in one run.
* Added `distmat` slot and `update` generic to save some time if possible. See examples of `dtwclust`.
* Extra parameters for distance functions should be correctly detected now.
* Using `dtw_lb` function now correctly warns about `pam.precompute` being `TRUE`.
* Option to calculate pairwise distances with `DTW_LB`, `SBD`, `LB_Keogh` and `LB_Improved`. See their respective notes. The distance function created for the `dtwclustFamily` slot also supports this.
* Clusters are randomly re-initialized if they become empty at some iteration.
* Now all included centroid functions recompute centers only if necessary.
* Option to optimize distmat calculation if the distance function is symmetric.

## Version 1.3.0
* Added the possibility to run several repetitions for partitional procedures, using different random starts each time by using the `doRNG` package
* Parallel computing all around. Watch out for RAM! Use Linux if possible.
     + Repetitions can be done in parallel if the user provides a suitable backend
     + Calculation of distance matrices also takes advantage of parallel backends.
     + `TADPole` and `dtw_lb` too, probably negligible (maybe even detrimental) for small datasets.

## Version 1.2.0
* Added the option to prevent pre-computation of distance matrix when using PAM centroids
* Added an example with a custom distance function
* Using closures instead of relying on passing environments as attributes
* Optimized the `SBD` function registered with `proxy`, it's a lot faster now
* Optimized clustering with `DBA` and `shape_extraction` so that centers are only recomputed if necessary
* Added processing time slot to class definition
* Several bug fixes (especially in case custom distances were used, and also for plot method)
* Fixed minor bugs in `TADPole`

## Version 1.1.0
* Added more options to the plot method for custom time labels
* Fixed an error in the calculation of Lemire's improved lower bound
* Added L2 norm to `DBA`, as well as window constraint
* Optimized `DBA` calculations
* More examples

## Version 1.0.0
* Limited support for time series of different lengths
* Improved plot method
* More control parameters
* Added more slots to the `dtwclust` class
* Some parameters changed order
* Default values of some auxiliary functions were changed (`DBA` and `SBD`)
* No longer supporting native `kccaFamilies`, supporting any distance registered in `proxy` instead
* Several bug fixes, especially if `dtw` or `dtw2` were being used, in which case they may have been erroneously computed before

## Version 0.1.0
* Initial release
