\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 4.0.0}{
    \itemize{
        \item Optimized \code{TADPole} for multiple \code{k} and \code{dc} values.
        \item Optimized PAM centroids with \code{pam.precompute = FALSE} by using sparse matrices from the \pkg{Matrix} package.
        \item Optimized \code{shape_extraction} by using the \code{eigs_sym} function from the \pkg{RSpectra} package.
        \item Implemented the DTW lower bounds in \code{C++}.
        \item Implemented DBA in \code{C++}.
        \item Implemented an alternative version of \emph{multivariate} DBA that might be faster. See its documentation.
        \item Added a \code{symmetric} control for fuzzy clustering.
        \item Partitional, hierarchical and fuzzy configurations in \code{compare_clusterings} now take into account the \code{symmetric} control parameter.
        \item The functionality for \code{pick.clus} in \code{compare_clusterings} changed depending on the value of \code{return.objects}.
        \item Fixed an error that sometimes caused objects returned by \code{tsclust} to have duplicated elements in the \code{args} slot.
        \item Fixed DTW symmetry detection for fuzzy clustering.
        \item Some internal functions changed, so older objects might no longer be compatible. Try using \code{update(old_TSClusters_obj)}.
        \item The \code{dtwclust} \emph{function} is now deprecated. Try using \code{as(dtwclust_class_obj, "TSClusters")} for old objects.
        \item Changed name of function \code{compute_envelop} to \code{compute_envelope} (old one still available but deprecated).
        \item Minor vignette updates.
        \item Several internal optimizations.
    }
}
