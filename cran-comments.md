## Update to version 2.0.0
* Refactored most of the code
* Dropped some formal parameters. Still supporting them via ellipsis.
* More parallel support and extended functionality. 
* New classes and some new generics.
* Bug fixes.

## Test environments
* Local Linux Mint 17.3, R 3.2.3
* Local Windows 10, R 3.2.3
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs, WARNINGs or NOTEs
* The testthat package has its own definition of 'proc_time' class, which results in several messages being displayed when the tests are run. This only happens during tests.
