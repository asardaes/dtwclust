# NEWS

## Version 1.1.1
* Added the option to prevent pre-computation of distance matrix when using PAM centroids
* Added an example with a custom distance function
* Using closures instead of relying on passing environments as attributes
* Bug fixes (in case custom distances were used)

## Version 1.1.0
* Added more options to the plot method for custom time labels
* Fixed an error in the calculation of Lemire's improved lower bound
* Added L2 norm to DBA, as well as window constraint
* Optimized DBA calculations
* More examples

## Version 1.0.0
* Limited support for time series of different lengths
* Improved plot method
* More control parameters
* Added more slots to the dtwclust class
* Some parameters changed order
* Default values of some auxiliary functions were changed (DBA and SBD)
* No longer supporting native kccaFamilies, supporting any distance registered in 'proxy' instead
* Several bug fixes, especially if 'dtw' or 'dtw2' were being used, in which case they may have been erroneously computed before

## Version 0.1.0
* Initial release
