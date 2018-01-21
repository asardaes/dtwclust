\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.1.0.9000}{
\itemize{
    \item Several functions now use \code{RcppParallel} for parallelization, refer to the new parallelization vignette for more information.
    \item None of the distance/centroid functions expose a parameter for helper matrices anymore.
    \item The ellipsis in \code{sdtw_cent} is now actually passed to \code{nloptr}.
    \item Fixed \code{tsclust} for \code{dtw_lb} distance and \code{pam.precompute = FALSE} (it was always using \code{dtw_basic}, so the results should not change).
    \item Updated distance and centroid timing experiments.
}
}
