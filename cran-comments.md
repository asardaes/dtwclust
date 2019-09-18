
## Update to version 5.5.5
- Fixed building of vignette accounting for new subcaption version.
- Added citation information, although the DOI page is not yet active.

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
