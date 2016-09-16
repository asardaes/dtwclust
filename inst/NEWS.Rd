\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.2.2}{
\itemize{
     \item Consistency adjustments: all \code{center(s)} arguments/slots will be removed and replaced with \code{centroid(s)}.
     \item The \code{SBD} and \code{shape_extraction} functions used the \code{crossprod} function. This has been changed, but it resulted in very small numerical differences, which may be enough to alter some clustering results.
}
}
