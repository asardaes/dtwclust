# Changelog

## Version 4.1.0
* Added some fuzzy cluster validity indices to `cvi`.
* When using `TADPole` clustering through `tsclust`, custom parameters can now be passed to a custom centroid function if provided.

## Version 4.0.3
* Deactivating PAM precomputation will no longer use sparse matrices by default, explicitly set `pam.sparse` to `TRUE` if you want this functionality.
* Fixed multivariate plots (#18).
* Fixed `zscore` for data frame input (it was still coercing column-wise instead of row-wise).
* Added an additional vignette with some timing experiments.

## Version 4.0.2
* Ported the `proxy` versions of `lb_keogh`, `lb_improved`, `SBD`, `GAK` and `dtw_basic` to `C++` and improved them by using the `bigmemory` package.
* Ported `TADPole` to `C++`.
* Ported part of the algorithm that updates sparse distance matrices to `C++`.
* Improved the optimizations for symmetric matrices that are calculated in parallel.
* Fixed `tsclustFamily`'s `dist` function for matrix or data frame input.
* Fixed partitional PAM centroids for `dtw_lb` distance and `pam.precompute = FALSE` (#16).
* Exported a function to coerce matrices and data frames to a list.
* Updated documentation.

## Version 4.0.1
* Ported `dtw_lb` to `C++` when using `dtw_basic`.
* Modified some tests to account for rounding error (CRAN request).

## Version 4.0.0
* Optimized `TADPole` for multiple `k` and `dc` values.
* Optimized PAM centroids with `pam.precompute = FALSE` by using sparse matrices from the `Matrix` package.
* Optimized `shape_extraction` by using the `eigs_sym` function from the `RSpectra` package.
* Implemented the DTW lower bounds in `C++`.
* Implemented DBA in `C++`.
  + Implemented an alternative version of *multivariate* DBA that might be faster. See its documentation.
* Added a `symmetric` control for fuzzy clustering.
* Partitional, hierarchical and fuzzy configurations in `compare_clusterings` now take into account the `symmetric` control parameter.
* The functionality for `pick.clus` in `compare_clusterings` changed depending on the value of `return.objects`.
* Fixed an error that sometimes caused objects returned by `tsclust` to have duplicated elements in the `args` slot.
* Fixed DTW symmetry detection for fuzzy clustering.
* Some internal functions changed, so older objects might no longer be compatible. Try using `update(old_TSClusters_obj)`.
* The `dtwclust` *function* is now deprecated. Try using `as(dtwclust_class_obj, "TSClusters")` for old objects.
* Changed name of function `compute_envelop` to `compute_envelope` (old one still available but deprecated).
* Minor vignette updates.
* Several internal optimizations.

## Version 3.2.0
* Added functions `compare_clusterings` and helpers to compare many clustering configurations, possibly in parallel.
* Fixed an error in `tsclust` that prevented CVIs to be calculated for hierarchical/TADPole cases if a custom centroid function was used.
* Added slot `seed` to the objects returned by `tsclust`.

## Version 3.1.2
* The arguments in `tsclust`'s ellipsis are now passed to all preprocessing, centroid and distance functions.
* Fixed some symmetry detection in `tsclust` when using DTW.
* Updated vignette to use `tsclust` in the examples.

## Version 3.1.1
* Seeds are now set when calling TADPole through `dtwclust` with a parallel backend and multiple values of `k`, in case the centroid function has randomness associated.
* Added a new experimental function `tsclust` that should be functionally equivalent to `dtwclust` but is hopefully more coherent in general.

## Version 3.1.0
* Fixed subsetting of multivariate plots
* Implemented the memory-saving version of DTW in `dtw_basic` when `backtrack = FALSE`.
* Implemented fuzzy c-medoids.
* Seeds were not being set in hierarchical/TADPole clustering, which could affect reproducibility of preprocessing/centroid functions.
* The `force.symmetry` argument in `dtw_lb` was removed since it didn't serve any real purpose.
* Exported the function to compute envelops.
* Minor vignette updates.

## Version 3.0.0
* Bear in mind that the `DTW/SBD` algorithms (and hence the functions that depend on them) might give different results in installations with 32-bit architectures.
* Removed deprecated arguments/slots.
    + If you have older `dtwclust` objects saved, try updating them with `attr(dtwclust_object, "centers") <- NULL` if you run into compatibility problems.
    + `DBA` arguments changed order.
* Removed the `dba.alignment` argument from `DBA` since other `step.pattern`s don't really work.
* Added (conditional) support for more hierarchical procedures. See the examples and vignette.
* Added a new distance based on global alignment kernels: `GAK`.
* Added support for functions in package `clue`.
* Fixed detection of some symmetric DTW cases.
* No longer enforcing preprocessing/centroid functions with ellipsis in their formals.
* Added a multivariate dataset sample: `CharTrajMV`.
* Improved plots for clusterings with multivariate series.
* Updated vignette.

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
