\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 4.0.0}{
    \itemize{
        \item Optimized \code{TADPole} for multiple \code{k} and \code{dc} values.
        \item Partitional and hierarchical configurations in \code{compare_clusterings} now take into account the \code{symmetric} control parameter if provided.
        \item Optimized PAM centroids with \code{pam.precompute = FALSE} by using sparse matrices from the \pkg{Matrix} package.
        \item Fixed an error that sometimes caused objects returned by \code{tsclust} to have duplicated elements in the \code{args} slot.
        \item Some internal functions changed, so older objects might no longer be compatible. Try using \code{update(old_obj)}.
    }
}
