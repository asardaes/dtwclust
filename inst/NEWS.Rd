\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.2.0.9000}{
\itemize{
    \item Added a shiny app, see \code{interactive_clustering}.
    \item Added another helper for \code{compare_clusterings}: \code{repeat_clustering}.
    \item Added \code{distmat} to \code{fuzzy_control} for when \code{centroid = "fcmdd"}.
    \item Dropped dependency on \code{rngtools}. Changing the random number generator after loading \code{dtwclust} should not be a problem now.
}
}
