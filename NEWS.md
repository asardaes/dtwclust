# NEWS

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
* Several hierarchical procedures can be made in one run.
* Added `distmat` slot and `update` generic to save some time if possible. See examples of `dtwclust`.
* Extra parameters for distance functions should be correctly detected now.
* Using `dtw_lb` function now correctly warns about `pam.precompute` being `TRUE`.
* Option to calculate pairwise distances with `SBD`, `LB_Keogh` and `LB_Improved`. See their respective notes. The distance function created for the `dtwclustFamily` slot also supports this.
* Clusters are randomly re-initialized if they become empty at some iteration.
* Now all included centroid functions recompute centers only if necessary.

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
