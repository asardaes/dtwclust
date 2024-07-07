
## Update to version 6.0.0
* Update Makevars for ARM version of Windows.
* Internal optimizations of do.call usages.

## Test environments
* Local GNU/Linux, R release
* win-builder (devel and release)
* GitHub CI
  + Linux: devel and release
  + OSX: release
  + Windows: release

## R CMD check results
* There were no ERRORs or WARNINGs
* There were 2 NOTEs:
  + The installed size is due to the 3 included vignettes (already compacted)
    and the compiled code (approx. 20% of the code).
  + GNU make is a system requirement due to RcppParallel
