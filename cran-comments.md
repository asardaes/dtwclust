
## Update to version 5.2.0
* Vignette updates, addition of multi-threading, and minor fixes.

## Test environments
* Local GNU/Linux, R release
* Local Windows 10, R release
* win-builder (devel and release)
* Travis CI
  + Linux: devel and release
  + OSX: release
* AppVeyor (x32 and x64)

## R CMD check results
* There were no ERRORs or WARNINGs
* There were 3 NOTEs:
  + Regarding spelling, TADPole is the name of an algorithm,
    and it is written as such. Centroid/Partitional are maybe domain specific, 
    but are written like that too.
  + The installed size is due to the 3 included vignettes (already compacted)
    and the compiled code (approx. 20% of the code).
  + GNU make is a system requirement due to RcppParallel

## Additional issues
* UBSAN: the warnings come from the bigmemory package.
  Apparently the sanitizers don't play nicely with its classes' inheritance.
* MKL: the equality tests have been adjusted.
