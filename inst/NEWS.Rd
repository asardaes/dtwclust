\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 6.0.0}{
\itemize{
    \item Update Makevars for ARM version of Windows.
    \item Sanitize internal usage of \code{do.call} to avoid huge backtraces.
    \item Support lower triangular \code{distmat} objects for symmetric distances. This is a slightly breaking change, see PR #77 on GitHub.
}
}
