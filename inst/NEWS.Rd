\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.0.0}{
    \itemize{
        \item Removed all deprecated code.
        \item Added external fuzzy cluster validity indices.
        \item Fixed GAK's \code{sigma} estimation in the non-proxy version. The estimation routine in the proxy version also changed slightly, and the estimates will be different due to randomness.
        \item Fixed the calculation of the Silhouette CVI.
        \item Modified \code{tslist} to hopefully make it faster.
        \item Partitional and fuzzy controls (\code{tsclust-controls}) now specify version 2 of the main loop by default.
        \item Packages \code{clue} and \code{ggplot2} are no longer automatically attached by default.
        \item Fixed minor potential bugs in several functions.
        \item Changed \code{plyr} dependency for \code{dplyr}.
        \item Minor vignette updates.
    }
}
