\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.3.0}{
\itemize{
     \item Exported the \code{dtw2} function.
     \item Fixed a bug where multivariate series with different length had spurious data added to them.
     \item Consistency adjustments: all \code{center(s)} arguments/slots will be removed and replaced with \code{centroid(s)}.
     \item Added the \code{dtw_basic} function, which might be marginally faster due to its limited functionality. It can also be used with \code{dtw_lb} and \code{DBA}. It will now be used by default, but in the case of \code{DBA}, the results \strong{will} be different because the backtracking algorithm differs with respect to the one in the \pkg{dtw} package.
     \item The \code{SBD} and \code{shape_extraction} functions used the \code{crossprod} function internally, which returned a 1x1 matrix by default, causing some dimension inconsistencies. In the future, \code{R} will give an error about the inconsistency, so the function has been changed in \pkg{dtwclust}, but it resulted in very small numerical differences, which may be enough to alter some clustering results.
}
}
