\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 3.0.0}{
    \itemize{
        \item Bear in mind that the \code{DTW/SBD} algorithms (and hence the functions that depend on them) might give different results in installations with 32-bit architectures.
        \item Removed deprecated arguments/slots. If you have older \code{dtwclust} objects saved, try updating them with \code{attr(dtwclust_object, "centers") <- NULL} if you run into compatibility problems.
        \item \code{DBA} arguments changed order.
        \item Removed the \code{dba.alignment} argument from \code{DBA} since other \code{step.pattern}s don't really work.
        \item Added (conditional) support for more hierarchical procedures. See the examples and vignette.
        \item Added a new distance based on global alignment kernels: \code{GAK}.
        \item Added support for functions in package \pkg{clue}.
        \item Fixed detection of some symmetric DTW cases.
        \item No longer enforcing preprocessing/centroid functions with ellipsis in their formals.
        \item Added a multivariate dataset sample: \code{CharTrajMV}.
        \item Improved plots for clusterings with multivariate series.
        \item Updated vignette.
    }
}
