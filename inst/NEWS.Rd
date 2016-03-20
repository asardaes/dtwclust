\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.1.1}{
\itemize{
     \item Minor bug fixes for fuzzy clustering
     \itemize{
          \item Fixed the \code{predict} generic.
          \item The final values returned in the \code{fcluster} slot needed one final update during clustering. It should be correct now, but it will vary slightly with respect to what was previously given in version 2.1.0.
     }
     \item Some more examples
}
}
