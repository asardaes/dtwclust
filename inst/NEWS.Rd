\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.2.1}{
\itemize{
     \item Added package vignette.
     \item Implemented several cluster validity indices in the new \code{cvi} function.
     \item The custom argument \code{force.pairwise} used in some of the \code{proxy} distances is not necessary anymore.
     \item No longer enforcing \code{NULL} dimensions for each series (in case multivariate series are to be used), use with caution.
     \item Fixed an error that prevented fuzzy clustering from working when \code{k = 2}.
     \item Fixed a possible problem in some of the included proxy distances when a data object had a length of 1.
}
}
