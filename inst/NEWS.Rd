\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 2.0.0}{
\itemize{
     \item Major refactor.
     \item Many formal parameters from the \code{dtwclust} function were dropped and implemented in the formal class \code{dtwclustControl}. For now, they will still be supported through \code{...} with a message.
     \item No longer supporting non-proxy distances, in order to be able to take advantage of the included optimizations.
     \item Dropped inheritance of \code{flexclust}'s \code{kccasimple}
     \itemize{
          \item Many slots and methods were ported
          \item Inheriting from \code{hclust} class now, and all its associated methods
     }
     \item \code{shape_extraction} now accepts series with different lengths!
     \item More parallel support
     \itemize{
          \item DBA and Shape centroid calculations
          \item DBA itself (probably unnecessary unless you're averaging a lot of series)
		\item Custom \code{proxy} distances directly (except \code{DTW2})
     }
     \item Several hierarchical procedures can be made in one run.
     \item Added \code{distmat} slot and \code{update} generic to save some time if possible. See examples of \code{dtwclust}
     \item Extra parameters for distance functions should be correctly detected now.
     \item Using \code{dtw_lb} function now correctly warns about \code{pam.precompute} being \code{TRUE}.
     \item Option to calculate pairwise distances with \code{DTW_LB}, \code{SBD}, \code{LB_Keogh} and \code{LB_Improved}. See their respective notes. The distance function created for the \code{dtwclustFamily} slot also supports this.
     \item Clusters are randomly re-initialized if they become empty at some iteration.
     \item Now all included centroid functions recompute centers only if necessary.
     \item Option to optimize distmat calculation if the distance function is symmetric.
}
}
