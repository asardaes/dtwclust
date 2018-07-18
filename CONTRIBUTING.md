# Contributing a distance measure

Thanks to the `proxy` package,
implementing a custom distance,
either entirely in R or using C/C++,
is relatively easy.
However, thanks to `RcppParallel`,
implementing a C++ distance with support for multi-threading is also possible.
The framework used by `dtwclust` will be explained here,
assuming familiarity with R's C/C++ interface and C++ itself.

## Foreword

One of the packages that is linked against is `RcppArmadillo`.
When one works with it,
the `Rcpp` header should *not* be included,
since its functionality is included automatically by `RcppArmadillo`.
However, the `RcppArmadillo` header increases compilation times considerably,
so it should only be included when necessary.
This is why many `dtwclust` source files contain several classes in the same file,
and it will also explain some of the design choices below;
this reduces compilation times significantly.

## Stand-alone function

Stand-alone means that it can be used directly without going through `proxy`.
On the R side there are pretty much no restrictions,
the function can have any number of parameters,
but it's probably better if consistency checks are done in R.
See for example [`dtw_basic`](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R).

On the C++ side the distance should be declared in the corresponding [header](https://github.com/asardaes/dtwclust/blob/master/src/distances/R-gateways.h),
registered in the [initialization](https://github.com/asardaes/dtwclust/blob/master/src/init.cpp),
and [defined](https://github.com/asardaes/dtwclust/blob/master/src/distances/R-gateways.cpp).
Importantly, if the stand-alone version will serve as a basis for the `proxy` version,
the core calculations should be done independently of any R/Rcpp API,
depending either on raw pointers,
or on some kind of wrapper;
such a wrapper is available in [this template](https://github.com/asardaes/dtwclust/blob/master/src/utils/SurrogateMatrix.h).
These core functions should be declared in the [internal distances header](https://github.com/asardaes/dtwclust/blob/master/src/distances/details.h).
Since they don't include many third-party headers,
they can be implemented in different files;
see for example [`dtw_basic`](https://github.com/asardaes/dtwclust/blob/master/src/distances/dtw-basic.cpp).

## `proxy` function

A [`DistanceCalculator`](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.h) has to be implemented.
The concrete class should be declared there,
and added to the [factory](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L28).
Since all time series are passed from R as a list of vectors, matrices, or complex vectors,
there is the `TSTSList` [templated class](https://github.com/asardaes/dtwclust/blob/master/src/utils/TSTSList.h) that works with `Rcpp`'s `NumericVector`, `NumericMatrix` and `ComplexVector`,
saving them as Armadillo's `mat` and `cx_mat` so that they are thread-safe.
Univariate series are saved as matrices with 1 column.

### Constructor

It is expected that the `DistanceCalculator`'s constructor will take 3 `SEXP` parameters that will all contain `Rcpp::List`s with:
any arguments for the distance,
the time series in `x` (from `proxy`),
and the time series in `y` (from `proxy`).
See for example the [DtwBasicCalculator](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L50),
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
For example, see what [DtwBasicCalculator](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L74) does.

### Calculate

The factory method `calculate` takes 2 integers `i` and `j`,
and it is expected that it returns the distance between `x[i]` and `y[j]`.
Some calculators [dispatch](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L69) to the appropriate method directly,
but [others](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L166) can pass more parameters as needed.
Also note how in most cases the core calculations are further [delegated](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L85) to the core functions described [above](#stand-alone-function),
but it's not always necessary.
For instance, the `SbdCalculator` does the [calculations directly](https://github.com/asardaes/dtwclust/blob/master/src/distances/calculators.cpp#L265),
since the stand-alone version is implemented entirely in R.

### R side

All functions have a similar structure (`dtw_basic` used as an example):

- [Check consistency of inputs](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L111)
- [Prepare common parameters by evaluationg a pre-defined expression](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L123)
- [Check consistency of distance parameters and put them in a list](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L126)
- [Get the available number of threads](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L155)
- [`.Call` `C_distmat_loop`](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L156)
- [Final output adjustments](https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R#L160)

They are registered with `proxy` during [loading](https://github.com/asardaes/dtwclust/blob/master/R/pkg.R#L77),
and unregistered during [unload](https://github.com/asardaes/dtwclust/blob/master/R/pkg.R#L154).

## Final details

The distance should be added to the [internal globals](https://github.com/asardaes/dtwclust/blob/master/R/UTILS-globals-internal.R#L9),
and there it can also specify if it supports series of different length and/or multivariate series.

There should be [unit tests](https://github.com/asardaes/dtwclust/blob/master/tests/testthat/unit/distances.R) for stand-alone functions,
as well as [integration tests](https://github.com/asardaes/dtwclust/blob/master/tests/testthat/integration/proxy.R) for `proxy` distances.
