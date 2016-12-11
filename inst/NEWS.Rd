\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 3.0.1}{
    \itemize{
        \item Fixed subsetting of multivariate plots.
        \item Implemented the memory-saving version of DTW in \code{dtw_basic} when \code{backtrack = FALSE}.
        \item The \code{force.symmetry} in \code{dtw_lb} is deprecated since it didn't serve any real purpose.
    }
}
