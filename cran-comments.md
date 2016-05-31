## Update to version 2.1.2
* Documentation fixes
* Added functionality
* Minor improvements

## Test environments
* Local GNU/Linux, R 3.3.0
* Local Windows 10, R 3.3.0
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs
* There was a NOTE regarding spelling. TADPole is the name of an algorithm, and it is written as such.
* The testthat package has its own definition of 'proc_time' class, which results in several messages being displayed when the tests are run. This only happens during tests.
