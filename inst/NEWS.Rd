\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 3.1.0}{
    \itemize{
        \item Fixed subsetting of multivariate plots.
        \item Implemented the memory-saving version of DTW in \code{dtw_basic} when \code{backtrack = FALSE}.
        \item Implemented fuzzy c-medoids.
        \item Seeds were not being set in hierarchical/TADPole clustering, which could affect reproducibility of preprocessing/centroid functions.
        \item The \code{force.symmetry} argument in \code{dtw_lb} was removed since it didn't serve any real purpose.
        \item Exported the function to compute envelops.
        \item Minor vignette updates.
    }
}
