\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 3.2.0}{
    \itemize{
        \item Added functions \code{compare_clusterings} and helpers to compare many clustering configurations, possibly in parallel.
        \item Fixed an error in \code{tsclust} that prevented CVIs to be calculated for hierarchical/TADPole cases if a custom centroid function was used.
        \item Added slot \code{seed} to the objects returned by \code{tsclust}.
    }
}
