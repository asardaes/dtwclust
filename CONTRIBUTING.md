# Contributing a distance measure

Thanks to the `proxy` package,
implementing a custom distance,
either entirely in R or using C/C++,
is relatively easy.
However, thanks to `RcppParallel`,
implementing a C++ distance with support for multi-threading is also possible.
The framework used by `dtwclust` will be explained here,
assuming familiarity with R's C/C++ interface and C++ itself.

## Stand-alone function

By stand-alone we mean that it can be used directly without going through `proxy`.
On the R side there are pretty much no restrictions,
the function can have any number of parameters,
but it's probably better if consistency checks are done in R.
See for example [`dtw_basic`](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L52).

On the C++ side the distance should be declared in the corresponding [header](https://github.com/asardaes/dtwclust/blob/master/src/distances/distances.h#L14),
registered in the [initialization](https://github.com/asardaes/dtwclust/blob/master/src/init.cpp#L8),
and [defined](https://github.com/asardaes/dtwclust/blob/master/src/distances/dtw-basic.cpp).
Importantly, if the stand-alone version will serve as a basis for the `proxy` version,
the [core calculations](https://github.com/asardaes/dtwclust/blob/master/src/distances/dtw-basic.cpp#L89) should be done independently of any R/Rcpp API,
depending either on raw pointers,
[custom wrapper classes](https://github.com/asardaes/dtwclust/blob/master/src/utils/SurrogateMatrix.h),
or `RcppParallel`'s wrappers.

## `proxy` function

For this case, we will start with the C++ side.
A [`DistanceCalculator`](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.h) has to be implemented.
The concrete class should be declared there,
and added to the [factory](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L27).
Since all time series are passed from R as a list of vectors, matrices, or complex vectors,
there is a [templated class](https://github.com/asardaes/dtwclust/blob/master/src/utils/TSTSList.h) that works with `Rcpp`'s `NumericVector`, `NumericMatrix` and `ComplexVector`,
saving them respectively as `RcppParallel`'s `RVector<double>`, `RMatrix<double>`, or Armadillo's `cx_vec` so that they are thread-safe.

### Constructor

It is expected that the `DistanceCalculator`'s constructor will take 3 `SEXP` parameters that will contain `Rcpp::List`s with:
any arguments for the distance,
the time series in `x` (from `proxy`),
and the time series in `y` (from `proxy`).
See for example the [DtwBasicCalculator](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L47),
and note that it handles univariate and multivariate series by expecting the R side to specify what was passed.

**Important**: if the concrete class has any private members that should be unique to each thread,
they should *not* be defined in the constructor.
See next part.

### Clone

All concrete implementations should have a clone method that returns a pointer to a new instance.
This method will be called from the different threads,
and each thread will `delete` the clone when it's done.
If there are private members that should be unique to each thread
(e.g. dynamically allocated memory),
they should be set-up during cloning.
See what [DtwBasicCalculator](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L101) does,
and note its [destructor](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L78).

### Calculate

The factory method `calculate` takes 2 integers `i` and `j`,
and it is expected that it returns the distance between `x[i]` and `y[j]`.
Some calculators [dispatch](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L85) to the appropriate (univariate/multivariate) method,
but [others](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L382) can do some more adjustments before dispatching.

### R side

All functions have a similar structure:

- [Check consistency of inputs](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L111)
- [Prepare common parameters by evaluationg a pre-defined expression](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L123)
- [Check consistency of distance parameters and put them in a list](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L126)
- [Get the available number of threads](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L156)
- [`.Call` `C_distmat_loop`](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L157)
- [Final output adjustments](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L161)

They are registered with `proxy` during [attachment](https://github.com/asardaes/dtwclust/blob/master/R/pkg.R#L76),
and unregistered during [unload](https://github.com/asardaes/dtwclust/blob/master/R/pkg.R#L149).

## Final details

The distance should be added to the [internal globals](https://github.com/asardaes/dtwclust/blob/master/R/UTILS-globals-internal.R#L9),
and there it can also specify if it supports series of different length and/or multivariate series.

There should be [unit tests](https://github.com/asardaes/dtwclust/blob/master/tests/testthat/unit/distances.R) for stand-alone functions,
as well as [integration tests](https://github.com/asardaes/dtwclust/blob/master/tests/testthat/integration/proxy.R) for `proxy` distances.
