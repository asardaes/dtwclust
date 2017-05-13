## Update to version 4.0.0
Several optimizations and bug fixes.
Some of the new functionality requires new dependencies,
and some previous dependencies were dropped.

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
