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

On the C++ side the distance should be declared in the corresponding [header](https://github.com/asardaes/dtwclust/blob/master/src/distances/distances.h),
registered in the [initialization](https://github.com/asardaes/dtwclust/blob/master/src/init.cpp#L8),
and [defined](https://github.com/asardaes/dtwclust/blob/master/src/distances/dtw-basic.cpp).
Importantly, if the stand-alone version will serve as a basis for the `proxy` version,
the [core calculations](https://github.com/asardaes/dtwclust/blob/master/src/distances/dtw-basic.cpp#L88) should be done independently of any R/Rcpp API,
depending either on raw pointers,
[custom wrapper classes](https://github.com/asardaes/dtwclust/blob/master/src/utils/SurrogateMatrix.h),
or `RcppParallel`'s wrappers.
These can be declared in the [internal distances header](https://github.com/asardaes/dtwclust/blob/master/src/distances/distances-details.h).

## `proxy` function

For this case, we will start with the C++ side.
A [`DistanceCalculator`](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.h) has to be implemented.
The concrete class should be declared there,
and added to the [factory](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L27).
Since all time series are passed from R as a list of vectors, matrices, or complex vectors,
there is the `TSTSList` [templated class](https://github.com/asardaes/dtwclust/blob/master/src/utils/TSTSList.h) that works with `Rcpp`'s `NumericVector`, `NumericMatrix` and `ComplexVector`,
saving them Armadillo's `mat` and `cx_mat` so that they are thread-safe.
Univariate series are saved as matrices with 1 column.

### Constructor

It is expected that the `DistanceCalculator`'s constructor will take 3 `SEXP` parameters that will contain `Rcpp::List`s with:
any arguments for the distance,
the time series in `x` (from `proxy`),
and the time series in `y` (from `proxy`).
See for example the [DtwBasicCalculator](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L53),
and note that `x` and `y` are simply given to the `TSTSList` template.

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
See what [DtwBasicCalculator](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L88) does,
and note its [destructor](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L70).

### Calculate

The factory method `calculate` takes 2 integers `i` and `j`,
and it is expected that it returns the distance between `x[i]` and `y[j]`.
Some calculators [dispatch](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L78) to the appropriate method directly,
but [others](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L324) can pass more parameters as needed.
Also note how in some cases the core calculations are further [delegated](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L100) by calling special wrappers defined in their own [header](https://github.com/asardaes/dtwclust/blob/master/src/distances/distances-details.h),
which don't use R or Rcpp's API.
This avoids duplicating code that is used in stand-alone functions (described [above](#stand-alone-function)),
but is not always necessary;
for instance, the SbdCalculator does the [calculations directly](https://github.com/asardaes/dtwclust/blob/master/src/distance-calculators/distance-calculators.cpp#L330),
since the stand-alone version is implemented entirely in R.

### R side

All functions have a similar structure:

- [Check consistency of inputs](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L111)
- [Prepare common parameters by evaluationg a pre-defined expression](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L123)
- [Check consistency of distance parameters and put them in a list](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L126)
- [Get the available number of threads](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L155)
- [`.Call` `C_distmat_loop`](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L156)
- [Final output adjustments](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L160)

They are registered with `proxy` during [loading](https://github.com/asardaes/dtwclust/blob/master/R/pkg.R#L77),
and unregistered during [unload](https://github.com/asardaes/dtwclust/blob/master/R/pkg.R#L152).

## Final details

The distance should be added to the [internal globals](https://github.com/asardaes/dtwclust/blob/master/R/UTILS-globals-internal.R#L9),
and there it can also specify if it supports series of different length and/or multivariate series.

There should be [unit tests](https://github.com/asardaes/dtwclust/blob/master/tests/testthat/unit/distances.R) for stand-alone functions,
as well as [integration tests](https://github.com/asardaes/dtwclust/blob/master/tests/testthat/integration/proxy.R) for `proxy` distances.
