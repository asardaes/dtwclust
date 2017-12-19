\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.1.0}{
\itemize{
    \item Added implementations of soft-DTW distance and centroid.
    \item Fixed the calculation of the \code{clusinfo} slot for \code{TSClusters} class for clusterings that ended with empty clusters.
    \item Fixed logic in partitional/fuzzy clustering, which was preventing the usage of paralle backends when only one value of k and one repetition were made.
}
}
