## Update to version 2.2.2
* Bug fixes
* Consistency adjustments

## Test environments
* Local GNU/Linux, R 3.3.1
* Local Windows 10, R 3.3.1
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs
* There was a NOTE regarding spelling: TADPole is the name of an algorithm, and it is written as such. Centroid/Partitional are maybe domain specific, but are written like that too.
* The "testthat" package has its own definition of 'proc_time' class, so said class is duplicated during tests, and only in that situation.
