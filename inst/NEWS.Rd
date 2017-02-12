\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\section{Changes in version 3.1.1}{
    \itemize{
        \item Seeds are now set when calling TADPole through \code{dtwclust} with a parallel backend and multiple values of \code{k}, in case the centroid function has randomness associated.
        \item Added a new experimental function \code{tsclust} that should be functionally equivalent to \code{dtwclust} but is hopefully more coherent in general.
    }
}
