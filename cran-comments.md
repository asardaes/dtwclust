## Update to version 3.0.0
* Removed deprecated arguments/slots
* Added support for more hierarchical methods
* Added support for functions in package clue
* Added a new distance measure

## Test environments
* Local GNU/Linux, R 3.3.2
* Local Windows 10, R 3.3.2
* win-builder (devel and release)
* Travis CI
  + Linux: devel and release
  + OSX: release
* AppVeyor

## R CMD check results
* There were no ERRORs or WARNINGs
* There was a NOTE regarding spelling: TADPole is the name of an algorithm, and it is written as such. Centroid/Partitional are maybe domain specific, but are written like that too.
