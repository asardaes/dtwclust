## Update to version 4.1.1
Short-term update: a modification introduced in v4.1.0 produces erroneous clustering results due to wrong logic.
Sorry for that.

## Test environments
* Local GNU/Linux, R 3.4.1
* Local Windows 10, R 3.4.1
* win-builder (devel and release)
* Travis CI
  + Linux: devel and release
  + OSX: release
* AppVeyor

## R CMD check results
* There were no ERRORs or WARNINGs
* There was a NOTE regarding spelling: TADPole is the name of an algorithm, and it is written as such. Centroid/Partitional are maybe domain specific, but are written like that too.
