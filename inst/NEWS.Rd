\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 4.1.1.9000}{
    \itemize{
        \item Removed all deprecated code.
        \item Fixed GAK's \code{sigma} estimation in the non-proxy version. The estimation routine in the proxy version also changed slightly.
        \item Modified \code{tslist} to hopefully make it faster.
        \item Partitional and fuzzy controls (\code{tsclust-controls}) now specify version 2 of the main loop by default.
        \item Packages \code{clue} and \code{ggplot2} are no longer automatically loaded by default.
        \item Fixed minor potential bugs in several functions.
    }
}
