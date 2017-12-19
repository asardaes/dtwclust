
## Update to version 5.1.0
* Adjusted tests for modified testthat API
* New distance and centroid function
* Minor bug fixes

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
* There was a NOTE regarding spelling: 
TADPole is the name of an algorithm, and it is written as such. 
Centroid/Partitional are maybe domain specific, but are written like that too.
