\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.3.0}{
\itemize{
    \item Added 2 shiny apps, see \code{interactive_clustering} and \code{ssdtwclust}.
    \item Added another helper for \code{compare_clusterings}: \code{repeat_clustering}.
    \item Added \code{distmat} to \code{fuzzy_control} for when \code{centroid = "fcmdd"}.
    \item Dropped dependency on \code{rngtools}. Changing the random number generator after loading \code{dtwclust} should not be a problem now.
    \item Fixed a bug regarding the way \code{tsclust} handled the \code{seed} when called from \code{compare_clusterings}. This will affect \code{compare_clusterings} results of hierarchical or TADPole clusterings made therein if a centroid subject to randomness was used (e.g. \code{dba}).
}
}
