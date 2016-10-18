\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 3.0.0}{
    \itemize{
        \item Removed deprecated arguments/slots. If you have older \code{dtwclust} objects saved, try updating them with \code{attr(dtwclust_object, "centers") <- NULL} if you run into compatibility problems.
        \item Added (conditional) support for more hierarchical procedures. See the examples.
        \item Added support for functions in package \pkg{clue}.
        \item Bear in mind that the \code{dtw_basic/SBD/shape_extraction} functions might give different results in 32-bit systems.
    }
}
