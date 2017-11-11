## Update to version 4.1.1.9000
* Removed deprecated code.
* New external fuzzy CVIs.
* Miscellaneous bug fixes and optimizations.

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
