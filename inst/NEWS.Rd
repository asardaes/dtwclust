\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.3.1.9000}{
\itemize{
    \item The \code{centroids} in the objects returned by \code{tsclust} now have an attribute \code{series_id} if said centroids were extracted from the provided series without modification. It is an integer vector indicating which series were chosen as centroids.
}
}
