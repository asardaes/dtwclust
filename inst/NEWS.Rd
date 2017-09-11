\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 4.1.1}{
    \itemize{
        \item Added some fuzzy cluster validity indices to \code{cvi}.
        \item When using \code{TADPole} clustering through \code{tsclust}, custom parameters can now be passed to a custom centroid function if provided.
        \item Modified the partitional/fuzzy loops, which were doing an extra iteration. Previous version still available for backwards compatibility. This modification was wrong in version 4.1.0.
    }
}
