\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.3.1.9000}{
\itemize{
    \item Fixed the calculation of SF and CH CVIs for hierarchical clustering with default centroid extraction.
    \item The \code{centroids} in the objects returned by \code{tsclust} now have an attribute \code{series_id} if said centroids were extracted from the provided series without modification. It is an integer vector indicating which series were chosen as centroids.
    \item Improved warning messages in \code{compare_clusterings}.
    \item Multi-threading parallelization is now more medium-grained to hopefully balance load a bit better.
}
}
