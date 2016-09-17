\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.2.2}{
\itemize{
     \item Consistency adjustments: all \code{center(s)} arguments/slots will be removed and replaced with \code{centroid(s)}.
     \item The \code{SBD} and \code{shape_extraction} functions used the \code{crossprod} function internally, which returned a 1x1 matrix by default, causing some dimension inconsistencies. In the future, \code{R} will give an error about the inconsistency, so the function has been changed in \pkg{dtwclust}, but it resulted in very small numerical differences, which may be enough to alter some clustering results.
}
}
