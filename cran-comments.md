## Update to version 4.0.1
* Short-term patch: modified some tests to account for rounding errors as requested (see https://www.r-project.org/nosvn/R.check/r-patched-solaris-x86/dtwclust-00check.html)

## Test environments
* Local GNU/Linux, R 3.4.0
* Local Windows 10, R 3.4.0
* win-builder (devel and release)
* Travis CI
  + Linux: devel and release
  + OSX: release
* AppVeyor

## R CMD check results
* There were no ERRORs or WARNINGs
* There was a NOTE regarding spelling: TADPole is the name of an algorithm, and it is written as such. Centroid/Partitional are maybe domain specific, but are written like that too.
