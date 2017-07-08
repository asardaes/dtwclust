\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 4.0.2}{
    \itemize{
        \item Ported the \code{proxy} versions of \code{lb_keogh}, \code{lb_improved}, \code{SBD}, \code{GAK} and \code{dtw_basic} to \code{C++} and improved them by using the \pkg{bigmemory} package.
        \item Ported \code{TADPole} to \code{C++}.
        \item Ported part of the algorithm that updates sparse distance matrices to \code{C++}.
        \item Improved the optimizations for symmetric matrices that are calculated in parallel.
        \item Fixed \code{tsclustFamily}'s \code{dist} function for matrix or data frame input.
        \item Fixed partitional PAM centroids for \code{dtw_lb} distance and \code{pam.precompute = FALSE} (issue #16 in GitHub).
        \item Exported a function to coerce matrices and data frames to a list.
        \item Updated documentation.
    }
}
