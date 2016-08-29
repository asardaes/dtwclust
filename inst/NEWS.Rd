\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.1.3}{
\itemize{
     \item Added package vignette.
     \item The custom argument \code{force.pairwise} used in some of the \code{proxy} distances is not necessary anymore.
     \item Added Silhouette Cluster Validity Index.
     \item No longer enforcing \code{NULL} dimensions for each series (in case multivariate series are to be used).
     \item Fixed an error that prevented fuzzy clustering to work when \code{k = 2}.
}
}
