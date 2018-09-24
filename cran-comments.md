
## Update to version 5.5.1
* Documentation fixes.
* Fixes related to setOldClass calls.

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
* There were 2 NOTEs:
  + The installed size is due to the 3 included vignettes (already compacted)
    and the compiled code (approx. 20% of the code).
  + GNU make is a system requirement due to RcppParallel

## Additional issues
* UBSAN: the errors stem from the RcppParallel package and can be safely ignored.
