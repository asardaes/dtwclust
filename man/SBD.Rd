% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DISTANCES-sbd.R
\name{SBD}
\alias{SBD}
\alias{sbd}
\title{Shape-based distance}
\usage{
SBD(x, y, znorm = FALSE, error.check = TRUE, return.shifted = TRUE)

sbd(x, y, znorm = FALSE, error.check = TRUE, return.shifted = TRUE)
}
\arguments{
\item{x, y}{Univariate time series.}

\item{znorm}{Logical. Should each series be z-normalized before calculating the distance?}

\item{error.check}{Logical indicating whether the function should try to detect inconsistencies
and give more informative errors messages. Also used internally to avoid repeating checks.}

\item{return.shifted}{Logical. Should the shifted version of \code{y} be returned? See details.}
}
\value{
For \code{return.shifted = FALSE}, the numeric distance value, otherwise a list with:
\itemize{
\item \code{dist}: The shape-based distance between \code{x} and \code{y}.
\item \code{yshift}: A shifted version of \code{y} so that it optimally matches \code{x} (based on \code{\link[=NCCc]{NCCc()}}).
}
}
\description{
Distance based on coefficient-normalized cross-correlation as proposed by Paparrizos and Gravano
(2015) for the k-Shape clustering algorithm.
}
\details{
This distance works best if the series are \emph{z-normalized}. If not, at least they should have
appropriate amplitudes, since the values of the signals \strong{do} affect the outcome.

If \code{x} and \code{y} do \strong{not} have the same length, it would be best if the longer sequence is
provided in \code{y}, because it will be shifted to match \code{x}. After matching, the series may have to
be truncated or extended and padded with zeros if needed.

The output values lie between 0 and 2, with 0 indicating perfect similarity.
}
\note{
If you wish to calculate the distance between several time series, it would be better to use the
version registered with the \code{proxy} package, since it includes some small optimizations. See the
examples.

This distance is calculated with help of the Fast Fourier Transform, so it can be sensitive to
numerical precision. Thus, this function (and the functions that depend on it) might return
different values in 32 bit installations compared to 64 bit ones.
}
\section{Proxy version}{


The version registered with \code{\link[proxy:dist]{proxy::dist()}} is custom (\code{loop = FALSE} in \link[proxy:registry]{proxy::pr_DB}).
The custom function handles multi-threaded parallelization directly with \link[RcppParallel:RcppParallel-package]{RcppParallel}.
It uses all available threads by default (see \code{\link[RcppParallel:setThreadOptions]{RcppParallel::defaultNumThreads()}}),
but this can be changed by the user with \code{\link[RcppParallel:setThreadOptions]{RcppParallel::setThreadOptions()}}.

An exception to the above is when it is called within a \code{\link[foreach:foreach]{foreach}} parallel loop \strong{made by dtwclust}.
If the parallel workers do not have the number of threads explicitly specified,
this function will default to 1 thread per worker.
See the parallelization vignette for more information - \code{browseVignettes("dtwclust")}

It also includes symmetric optimizations to calculate only half a distance matrix when appropriate---only one list of series should be provided in \code{x}.
Starting with version 6.0.0, this optimization means that the function returns an array with the lower triangular values of the distance matrix,
similar to what \code{\link[stats:dist]{stats::dist()}} does;
see \linkS4class{DistmatLowerTriangular} for a helper to access elements as it if were a normal matrix.
If you want to avoid this optimization, call \link[proxy:dist]{proxy::dist} by giving the same list of series in both \code{x} and \code{y}.

In some situations, e.g. for relatively small distance matrices, the overhead introduced by the
logic that computes only half the distance matrix can be bigger than just calculating the whole
matrix.
}

\examples{

# load data
data(uciCT)

# distance between series of different lengths
sbd <- SBD(CharTraj[[1]], CharTraj[[100]], znorm = TRUE)$dist

# cross-distance matrix for series subset (notice the two-list input)
sbD <- proxy::dist(CharTraj[1:10], CharTraj[1:10], method = "SBD", znorm = TRUE)

}
\references{
Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
In \emph{Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data}, series
SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \doi{10.1145/2723372.2737793}.
}
\seealso{
\code{\link[=NCCc]{NCCc()}}, \code{\link[=shape_extraction]{shape_extraction()}}
}
