#' Time series clustering with DTW
#'
#' This function uses the DTW distance and related lower bounds to cluster time series. For now, all series
#' must have equal length.
#'
#' @param data Numeric matrix where each row is a time series.
#' @param type What type of clustering method to use, \code{partitional}, \code{hierarchical} or \code{tadpole}.
#' @param k Numer of desired clusters in partitional methods.
#' @param method Which linkage method to use in hierarchical methods.
#' @param distance One of the supported distance measurements (see details). It can also be the name of a
#' family to use with function \code{\link[flexclust]{kcca}} if \code{type == "partitional"}, or a supported
#' distance of \code{\link[proxy]{dist}} if \code{type == "hierarchical"}.
#' @param centroid Either a supported string (see details) or an appropriate function to calculate centroids
#' when using partitional methods.
#' @param window.size Window constraint for DTW and LB calculations.
#' @param norm Pointwise distance. Either \code{L1} for Manhattan distance or \code{L2} for Euclidean.
#' @param dc Cutoff distance for TADPole algorithm.
#' @param seed Random seed for reproducibility
#' @param control Parameters for partitional clustering algorithm. See \code{\link[flexclust]{flexclustControl}}.
#' @param save.data Return a copy of the data in the returned object?
#' @param trace Boolean flag. If true, more output regarding the progress is printed to screen.
#' @param ... Additional arguments to pass to \code{\link[proxy]{dist}}.
#'
#' @return An object with formal class \code{\link{dtwclust-class}} if \code{type == "partitional" | "tadpole"}. Otherwise
#' an object with class \code{hclust} as returned by \code{\link[stats]{hclust}}.
#'
#' @author Alexis Sarda
#'
#' @export
#' @import flexclust
#' @importFrom proxy dist
#' @importFrom modeltools ModelEnvMatrix

dtwclust <- function(data = NULL, type = "partitional", k = 2, method = "average",
                     distance = "dtw_lb", centroid = "median",
                     window.size = NULL, norm = "L1", dc,
                     control = NULL, save.data = FALSE,
                     seed = NULL, trace = F,
                     ...)
{
     if (is.null(data))
          stop("No data provided")

     tic <- proc.time()

     type <- match.arg(type, c("partitional", "hierarchical", "tadpole"))
     norm <- match.arg(norm, c("L1", "L2"))

#      distance <- match.arg(distance, c("dtw", "dtw_lb", "lbk", "lbi", "sbd",
#                                        "kmeans", "kmedians", "angle", "jaccard", "ejaccard"))

     if (distance == "dtw") {
          distance <- switch(EXPR = norm,
                             L1 = "dtw",
                             L2 = "dtw2")
     }

     if (type == "partitional") {

          ## =================================================================================================================
          ## Partitional
          ## =================================================================================================================

          if (is.function(centroid))
               cent <- centroid
          else {
               centroid <- match.arg(centroid, c("mean", "median", "shape", "dba", "pam"))

               cent <- switch(EXPR = centroid,

                              mean = function(x) { return(apply(x, 2, mean)) },

                              median = function(x) { return(apply(x, 2, median)) },

                              shape = shape_extraction, # placeholder, will be ignored (see utils.R)

                              dba = DBA, # placeholder, will be ignored (see utils.R)

                              pam = function(x) NULL, # placeholder, will be ignored (see utils.R)

                              ## Otherwise
                              stop("Unsupported centroid calculation method"))
          }

          family <- switch(EXPR = distance,

                           ## Full DTW with L1 norm
                           dtw = kccaFamily(name = "DTW",
                                            dist = function(x, centers) {
                                                 if (!is.null(attr(x, "distmat")))
                                                      d <- dsub_pam(x, centers)
                                                 else if (is.null(window.size)) {
                                                      d <- proxy::dist(x = x, y = centers,
                                                                       method = "DTW", dist.method = "L1",
                                                                       ...)
                                                 } else {
                                                      d <- proxy::dist(x = x, y = centers,
                                                                       method = "DTW", dist.method = "L1",
                                                                       window.type = "sakoechiba",
                                                                       window.size = window.size,
                                                                       ...)
                                                 }

                                                 return(d)
                                            },


                                            cent = cent),

                           ## Full DTW with L2 norm
                           dtw2 = kccaFamily(name = "DTW2",
                                            dist = function(x, centers) {
                                                 if (!is.null(attr(x, "distmat")))
                                                      d <- dsub_pam(x, centers)
                                                 else if (is.null(window.size)) {
                                                      d <- proxy::dist(x = x, y = centers,
                                                                       method = "DTW2",
                                                                       ...)
                                                 } else {
                                                      d <- proxy::dist(x = x, y = centers,
                                                                       method = "DTW2",
                                                                       window.type = "sakoechiba",
                                                                       window.size = window.size,
                                                                       ...)
                                                 }

                                                 return(d)
                                            },


                                            cent = cent),


                           ## DTW with aid of lower bounds
                           dtw_lb = {
                                if (is.null(window.size))
                                     stop("You must provide a window size for this method")

                                kccaFamily(name = "DTW_LB",

                                           dist = function(x, centers) {
                                                if (!is.null(attr(x, "distmat")))
                                                     d <- dsub_pam(x, centers)
                                                else
                                                     d <- dtw_lb(x, centers, window.size, error.check=FALSE)

                                                return(d)
                                           },

                                           cent = cent)
                           },



                           ## Lemire's improved lower bound with L1
                           lbi = {
                                if (is.null(window.size))
                                     stop("You must provide a window size for this method")

                                kccaFamily(name = "LB_Improved",
                                           dist = function(x, centers) {
                                                if (!is.null(attr(x, "distmat")))
                                                     d <- dsub_pam(x, centers)
                                                else {
                                                     d <- proxy::dist(x = x, y = centers,
                                                                      method = "LBI",
                                                                      norm = norm,
                                                                      window.size = window.size,
                                                                      error.check=FALSE,
                                                                      ...)
                                                }

                                                return(d)
                                           },


                                           cent = cent)
                           },

                           ## Keogh's lower bound but with L1
                           lbk = {
                                if (is.null(window.size))
                                     stop("You must provide a window size for this method")

                                kccaFamily(name = "LB_Keogh",
                                           dist = function(x, centers) {
                                                if (!is.null(attr(x, "distmat")))
                                                     d <- dsub_pam(x, centers)
                                                else {
                                                     d <- proxy::dist(x = x, y = centers,
                                                                      method = "LBK",
                                                                      norm = norm,
                                                                      window.size = window.size,
                                                                      error.check=FALSE,
                                                                      ...)
                                                }

                                                return(d)
                                           },


                                           cent = cent)
                           },

                           ## Paparrizos' shape-based distance
                           sbd = {
                                kccaFamily(name = "SBD",
                                           dist = function(x, centers) {
                                                if (!is.null(attr(x, "distmat")))
                                                     d <- dsub_pam(x, centers)
                                                else {
                                                     d <- proxy::dist(x = x, y = centers,
                                                                      method = "SBD",
                                                                      ...)
                                                }

                                                return(d)
                                           },


                                           cent = cent)
                           },


                           ## Otherwise
                           kccaFamily(distance)
          )

          if (centroid == "shape") {
               family@allcent <- allcent_se

               family@preproc <- preproc_se

          } else if (centroid == "dba") {
               family@allcent <- allcent_dba

          } else if (centroid == "pam") {
               family@allcent <- allcent_pam

               family@preproc <- preproc_pam
          }


          if (is.null(control))
               ctrl <- new("flexclustControl")
          else
               ctrl <- as(control, "flexclustControl")

          if (trace)
               ctrl@verbose <- 1

          if (!is.null(seed))
               set.seed(seed)

          kc <- kcca(x = data,
                     k = k,
                     family = family,
                     simple = TRUE,
                     control = ctrl,
                     save.data = save.data)

          toc <- proc.time() - tic

          if (trace)
               cat("\n\tElapsed time is", toc["elapsed"], "seconds.\n\n")

          dtwc <- new("dtwclust", kc,
                      type = type,
                      distance = distance,
                      centroid = ifelse(is.function(centroid), "custom", centroid))

          dtwc

     } else if (type == "hierarchical") {

          ## =================================================================================================================
          ## Hierarchical
          ## =================================================================================================================

          if (trace)
               cat("\n\tCalculating distance matrix...\n")

          if (class(data) == "matrix")
               x <- lapply(seq_len(nrow(data)), function(i) data[i,])

          D <- switch(EXPR = distance,

                      ## Full DTW with L1 norm
                      dtw = {
                           if (is.null(window.size)) {
                                d <- proxy::dist(x = x, y = x,
                                                 method = "DTW", dist.method = "L1",
                                                 ...)
                           } else {
                                d <- proxy::dist(x = x, y = x,
                                                 method = "DTW", dist.method = "L1",
                                                 window.type = "sakoechiba", window.size = window.size,
                                                 ...)
                           }

                           d
                      },

                      ## Full DTW with L2 norm
                      dtw2 = {
                           if (is.null(window.size)) {
                                d <- proxy::dist(x = x, y = x,
                                                 method = "DTW2",
                                                 ...)
                           } else {
                                d <- proxy::dist(x = x, y = x,
                                                 method = "DTW2",
                                                 window.type = "sakoechiba", window.size = window.size,
                                                 ...)
                           }

                           d
                      },

                      ## DTW with aid of lower bounds
                      dtw_lb = {
                           if (is.null(window.size))
                                stop("You must provide a window size for this method")

                           d <- dtw_lb(x, x, window.size, error.check=TRUE)

                           d
                      },


                      ## Lemire's improved lower bound with L1
                      lbi = {
                           if (is.null(window.size))
                                stop("You must provide a window size for this method")

                           d <- proxy::dist(x = x, y = x,
                                            method = "LBI", window.size = window.size, norm = norm, error.check=TRUE,
                                            ...)

                           d
                      },

                      ## Keogh's lower bound but with L1
                      lbk = {
                           if (is.null(window.size))
                                stop("You must provide a window size for this method")

                           d <- proxy::dist(x = x, y = x,
                                            method = "LBK", window.size = window.size, norm = norm, error.check=TRUE,
                                            ...)

                           d
                      },

                      ## Paparrizos' shape-based distance
                      sbd = {
                           d <- proxy::dist(x = x, y = x,
                                            method = "SBD",
                                            ...)
                           d
                      },


                      ## Otherwise
                      proxy::dist(x = x, y = x, method = distance, ...)
          )

          if (trace)
               cat("\n\tPerforming hierarchical clustering...\n")

          ## Required form for 'hclust'
          D <- D[lower.tri(D)]

          ## Needed attribute for 'hclust' (case sensitive)
          attr(D, "Size") <- length(x)

          hc <- hclust(D, method = method)

          toc <- proc.time() - tic

          if (trace)
               cat("\n\tElapsed time is", toc["elapsed"], "seconds.\n\n")

          hc

     } else if (type == "tadpole") {

          ## =================================================================================================================
          ## TADPole
          ## =================================================================================================================

          if (is.null(window.size)) {
               stop("Please provide the 'window.size' parameter")
          }
          if (window.size%%2 != 0) {
               stop("For the Sakoe-Chiba band, the window must be symmetric and window.size must be even")
          }
          if (window.size <= 1) {
               stop("Window width must be larger than 1")
          }

          if (class(data) == "matrix")
               x <- lapply(seq_len(nrow(data)), function(i) data[i,])

          consistency_check(x, "tslist")

          R <- TADPole(x, window.size = window.size, k = k, dc = dc, error.check = FALSE)

          if (save.data) {
               tadpc <- new("dtwclust",
                            type = type,
                            distance = "dtw",
                            centroid = "TADPole (PAM)",

                            centers = data[R$centers, ],
                            k = as.integer(k),
                            cluster = as.integer(R$cl),
                            data = ModelEnvMatrix(designMatrix = data))
          } else {
               tadpc <- new("dtwclust",
                            type = type,
                            distance = "dtw",
                            centroid = "TADPole (PAM)",

                            centers = data[R$centers, ],
                            k = as.integer(k),
                            cluster = as.integer(R$cl))
          }

          toc <- proc.time() - tic

          if (trace)
               cat("\n\tElapsed time is", toc["elapsed"], "seconds.\n\n")

          tadpc
     }
}
