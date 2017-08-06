\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 4.0.3}{
    \itemize{
        \item Deactivating PAM precomputation will no longer use sparse matrices by default, explicitly set \code{pam.sparse} to \code{TRUE} if you want this functionality.
        \item Fixed multivariate plots.
        \item Fixed \code{zscore} for data frame input (it was still coercing column-wise instead of row-wise).
        \item Added an additional vignette with some timing experiments.
    }
}
