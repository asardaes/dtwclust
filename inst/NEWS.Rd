\name{NEWS}
\alias{NEWS}

\title{News for package \pkg{dtwclust}}

\description{
    Full changelog available at \url{https://github.com/asardaes/dtwclust/blob/master/CHANGELOG.md}
}

\section{Changes in version 5.5.7}{
\itemize{
    \item Series order is now preserved in plots of included shiny apps.
    \item The optimization in \code{sdtw_cent} now uses \code{stats::optim} by default.
    \item Dropped strong dependency on \code{bigmemory} in case it is archived. It can still be used if available.
}
}
