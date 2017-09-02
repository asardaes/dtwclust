#' Methods for `dtwclust`
#'
#' **Deprecated** methods associated with [dtwclust-class] objects.
#'
#' @name dtwclust-methods
#' @rdname dtwclust-methods
#' @include dtwclust-classes.R
#' @keywords internal
#'
#' @details
#'
#' Please refer to [tsclusters-methods] for the updated versions.
#'
NULL

# ==================================================================================================
# Custom initialize
# ==================================================================================================

## to avoid infinite recursion (see https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16629)
setMethod("initialize", "dtwclust",
          function(.Object, ..., call) {
              .Object <- methods::callNextMethod(.Object = .Object, ...)

              if (!missing(call)) .Object@call <- call

              .Object
          })

## for dtwclustFamily
setMethod("initialize", "dtwclustFamily",
          function(.Object, dist, allcent, ..., distmat = NULL, control = NULL, fuzzy = FALSE) {
              dots <- list(...)
              dots$.Object <- .Object

              control <- methods::as(control, "dtwclustControl")

              if (!missing(dist)) {
                  if (is.character(dist))
                      dots$dist <- ddist(dist, control, distmat)
                  else
                      dots$dist <- dist
              }

              if (fuzzy) {
                  dots$cluster <- fcm_cluster # fuzzy.R
                  if (!missing(allcent)) allcent <- match.arg(allcent, c("fcm", "fcmdd"))
              }

              if (!missing(allcent)) {
                  if (is.character(allcent))
                      dots$allcent <- all_cent(allcent,
                                               distmat = distmat,
                                               distfun = dots$dist,
                                               control = control,
                                               fuzzy = fuzzy)
                  else if (is.function(allcent))
                      dots$allcent <- allcent
                  else
                      stop("Centroid definition must be either a function or a character")
              }

              do.call(callNextMethod, dots)
          })

# ==================================================================================================
# Show
# ==================================================================================================

#' @rdname dtwclust-methods
#' @aliases show,dtwclust
#' @exportMethod show
#'
#' @param object,x An object of class [dtwclust-class] as returned by [dtwclust()].
#'
setMethod("show", "dtwclust",
          function(object) {
              print(object@call)
              cat("\n")

              cat(object@type, "clustering with", object@k, "clusters\n")
              cat("Using", object@distance, "distance\n")
              cat("Using", object@centroid, "centroids\n")

              if (object@type == "hierarchical")
                  cat("Using method", object@method, "\n")
              if (object@preproc != "none")
                  cat("Using", object@preproc, "preprocessing\n")

              cat("\nTime required for analysis:\n")
              print(object@proctime)

              if (object@type == "fuzzy") {
                  cat("\nHead of fuzzy memberships:\n\n")
                  print(utils::head(object@fcluster))

              } else {
                  cat("\nCluster sizes with average intra-cluster distance:\n\n")
                  print(object@clusinfo)
              }

              invisible(NULL)
          })

# ==================================================================================================
# update from stats
# ==================================================================================================

#' @rdname dtwclust-methods
#' @method update dtwclust
#' @export
#'
#' @param evaluate Logical. Defaults to `TRUE` and evaluates the updated call, which will result in
#'   a new `dtwclust` object. Otherwise, it returns the unevaluated call.
#'
update.dtwclust <- function(object, ..., evaluate = TRUE) {
    args <- as.pairlist(list(...))

    if (length(args) == 0L) {
        message("Nothing to be updated")

        if (evaluate)
            return(object)
        else
            return(object@call)
    }

    new_call <- object@call
    new_call[names(args)] <- args

    if (evaluate)
        ret <- eval.parent(new_call, n = 2L)
    else
        ret <- new_call

    ret
}

#' @rdname dtwclust-methods
#' @aliases update,dtwclust
#' @exportMethod update
#'
setMethod("update", methods::signature(object = "dtwclust"), update.dtwclust)

# ==================================================================================================
# predict from stats
# ==================================================================================================

#' @rdname dtwclust-methods
#' @method predict dtwclust
#' @export
#'
#' @param newdata New data to be assigned to a cluster. It can take any of the supported formats of
#'   [dtwclust()]. Note that for multivariate series, this means that it **must** be a list of
#'   matrices, even if the list has only one element.
#'
predict.dtwclust <- function(object, newdata = NULL, ...) {
    if (is.null(newdata)) {
        if (object@type != "fuzzy")
            ret <- object@cluster
        else
            ret <- object@fcluster

    } else {
        newdata <- tslist(newdata)
        check_consistency(newdata, "vltslist")
        nm <- names(newdata)

        newdata <- do.call(object@family@preproc,
                           args = enlist(newdata,
                                         dots = object@dots))

        distmat <- do.call(object@family@dist,
                           args = enlist(x = newdata,
                                         centroids = object@centroids,
                                         dots = object@dots))

        ret <- object@family@cluster(distmat = distmat, m = object@control@fuzziness)

        if (object@type != "fuzzy")
            names(ret) <- nm
        else
            dimnames(ret) <- list(nm, paste0("cluster_", 1L:ncol(ret)))
    }

    ret
}

#' @rdname dtwclust-methods
#' @aliases predict,dtwclust
#' @exportMethod predict
#'
setMethod("predict", methods::signature(object = "dtwclust"), predict.dtwclust)

# ==================================================================================================
# Plot
# ==================================================================================================

#' @rdname dtwclust-methods
#' @method plot dtwclust
#' @export
#'
#' @param y Ignored.
#' @param ... For `plot`, further arguments to pass to [ggplot2::geom_line()] for the plotting of
#'   the *cluster centroids*, or to [stats::plot.hclust()]. See details. For `update`, any supported
#'   argument. Otherwise ignored.
#' @param clus A numeric vector indicating which clusters to plot.
#' @param labs.arg Arguments to change the title and/or axis labels. See [ggplot2::labs()] for more
#'   information
#' @param data Optionally, the data in the same format as it was provided to [dtwclust()].
#' @param time Optional values for the time axis. If series have different lengths, provide the time
#'   values of the longest series.
#' @param plot Logical flag. You can set this to `FALSE` in case you want to save the ggplot object
#'   without printing anything to screen
#' @param type What to plot. `NULL` means default. See details.
#'
plot.dtwclust <- function(x, y, ...,
                          clus = seq_len(x@k), labs.arg = NULL,
                          data = NULL, time = NULL,
                          plot = TRUE, type = NULL)
{
    ## set default type if none was provided
    if (!is.null(type))
        type <- match.arg(type, c("dendrogram", "series", "centroids", "sc"))
    else if (x@type == "hierarchical")
        type <- "dendrogram"
    else
        type <- "sc"

    ## plot dendrogram?
    if (x@type == "hierarchical" && type == "dendrogram") {
        x <- S3Part(x, strictS3 = TRUE)
        if (plot) graphics::plot(x, ...)
        return(invisible(NULL))

    } else if (x@type != "hierarchical" && type == "dendrogram") {
        stop("Dendrogram plot only applies to hierarchical clustering.")
    }

    ## Obtain data, the priority is: provided data > included data list
    if (!is.null(data)) {
        data <- tslist(data)

    } else {
        if (length(x@datalist) < 1L)
            stop("Provided object has no data. Please re-run the algorithm with save.data = TRUE ",
                 "or provide the data manually.")

        data <- x@datalist
    }

    ## centroids consistency
    check_consistency(centroids <- x@centroids, "vltslist")

    ## force same length for all multivariate series/centroids in the same cluster by
    ## adding NAs
    if (mv <- is_multivariate(data)) {
        clusters <- split(data, factor(x@cluster, levels = 1L:x@k), drop = FALSE)

        for (id_clus in 1L:x@k) {
            cluster <- clusters[[id_clus]]
            if (length(cluster) < 1L) next ## empty cluster

            nc <- NCOL(cluster[[1L]])
            len <- sapply(cluster, NROW)
            L <- max(len, NROW(centroids[[id_clus]]))
            trail <- L - len

            clusters[[id_clus]] <- mapply(cluster, trail,
                                          SIMPLIFY = FALSE,
                                          FUN = function(mvs, trail) {
                                              rbind(mvs, matrix(NA, trail, nc))
                                          })

            trail <- L - NROW(centroids[[id_clus]])

            centroids[[id_clus]] <- rbind(centroids[[id_clus]], matrix(NA, trail, nc))
        }

        ## split returns the result in order of the factor levels,
        ## but I want to keep the original order as returned from clustering
        ido <- sort(sort(x@cluster, index.return=T)$ix, index.return = TRUE)$ix
        data <- unlist(clusters, recursive = FALSE)[ido]
    }

    ## helper values
    L1 <- lengths(data)
    L2 <- lengths(centroids)

    ## timestamp consistency
    if (!is.null(time) && length(time) < max(L1, L2))
        stop("Length mismatch between values and timestamps")

    ## Check if data was z-normalized
    if (x@preproc == "zscore")
        title_str <- "Clusters' members (z-normalized)"
    else
        title_str <- "Clusters' members"

    ## transform to data frames
    dfm <- reshape2::melt(data)
    dfcm <- reshape2::melt(centroids)

    ## time, cluster and colour indices
    color_ids <- integer(x@k)

    dfm_tcc <- mapply(x@cluster, L1, USE.NAMES = FALSE, SIMPLIFY = FALSE,
                      FUN = function(clus, len) {
                          t <- if (is.null(time)) seq_len(len) else time[1L:len]
                          cl <- rep(clus, len)
                          color <- rep(color_ids[clus], len)

                          color_ids[clus] <<- color_ids[clus] + 1L

                          data.frame(t = t, cl = cl, color = color)
                      })

    dfcm_tc <- mapply(1L:x@k, L2, USE.NAMES = FALSE, SIMPLIFY = FALSE,
                      FUN = function(clus, len) {
                          t <- if (is.null(time)) seq_len(len) else time[1L:len]
                          cl <- rep(clus, len)

                          data.frame(t = t, cl = cl)
                      })

    ## bind
    dfm <- data.frame(dfm, do.call(rbind, dfm_tcc))
    dfcm <- data.frame(dfcm, do.call(rbind, dfcm_tc))

    ## make factor
    dfm$cl <- factor(dfm$cl)
    dfcm$cl <- factor(dfcm$cl)
    dfm$color <- factor(dfm$color)

    ## create gg object
    gg <- ggplot2::ggplot(data.frame(t = integer(),
                                     variable = factor(),
                                     value = numeric(),
                                     cl = factor(),
                                     color = factor()),
                          aes_string(x = "t",
                                     y = "value",
                                     group = "L1"))

    ## add centroids first if appropriate, so that they are at the very back
    if (type %in% c("sc", "centroids")) {
        if (length(list(...)) == 0L)
            gg <- gg + ggplot2::geom_line(data = dfcm[dfcm$cl %in% clus, ],
                                          linetype = "dashed",
                                          size = 1.5,
                                          colour = "black",
                                          alpha = 0.5)
        else
            gg <- gg + ggplot2::geom_line(data = dfcm[dfcm$cl %in% clus, ], ...)
    }

    ## add series next if appropriate
    if (type %in% c("sc", "series")) {
        gg <- gg + ggplot2::geom_line(data = dfm[dfm$cl %in% clus, ],
                                      aes_string(colour = "color"))
    }

    ## add vertical lines to separate variables of multivariate series
    if (mv) {
        ggdata <- data.frame(cl = rep(1L:x@k, each = (nc - 1L)),
                             vbreaks = as.numeric(1L:(nc - 1L) %o% sapply(centroids, NROW)))

        gg <- gg + ggplot2::geom_vline(data = ggdata[ggdata$cl %in% clus, , drop = FALSE],
                                       colour = "black", linetype = "longdash",
                                       aes_string(xintercept = "vbreaks"))
    }

    ## add facets, remove legend, apply kinda black-white theme
    gg <- gg +
        ggplot2::facet_wrap(~cl, scales = "free_y") +
        ggplot2::guides(colour = FALSE) +
        ggplot2::theme_bw()

    ## labels
    if (!is.null(labs.arg))
        gg <- gg + ggplot2::labs(labs.arg)
    else
        gg <- gg + ggplot2::labs(title = title_str)

    ## plot without warnings in case I added NAs for multivariate cases
    if (plot) suppressWarnings(graphics::plot(gg))

    invisible(gg)
}

#' @rdname dtwclust-methods
#' @aliases plot,dtwclust,missing
#' @exportMethod plot
#'
setMethod("plot", methods::signature(x = "dtwclust", y = "missing"), plot.dtwclust)

# ==================================================================================================
# Cluster validity indices
# ==================================================================================================

#' @rdname cvi
#' @aliases cvi,dtwclust
#' @exportMethod cvi
#'
setMethod("cvi", methods::signature(a = "dtwclust"),
          function(a, b = NULL, type = "valid", ...) {
              if (a@type == "fuzzy")
                  stop("Only CVIs for crisp partitions are currently implemented.")

              type <- match.arg(type, several.ok = TRUE,
                                c("RI", "ARI", "J", "FM", "VI",
                                  "Sil", "SF", "CH", "DB", "DBstar", "D", "COP",
                                  "valid", "internal", "external"))

              dots <- list(...)

              internal <- c("Sil", "SF", "CH", "DB", "DBstar", "D", "COP")
              external <- c("RI", "ARI", "J", "FM", "VI")

              if (any(type == "valid")) {
                  type <- if (is.null(b)) internal else c(internal, external)

              } else if (any(type == "internal")) {
                  type <- internal

              } else if (any(type == "external")) {
                  type <- external
              }

              which_internal <- type %in% internal
              which_external <- type %in% external

              if (any(which_external))
                  CVIs <- cvi(a@cluster, b = b, type = type[which_external], ...)
              else
                  CVIs <- numeric()

              type <- type[which_internal]

              if (any(which_internal)) {
                  if (length(a@datalist) == 0L && any(type %in% c("SF", "CH"))) {
                      warning("Internal CVIs: the original data must be in object to calculate ",
                              "the following indices:",
                              "\n\tSF\tCH")

                      type <- setdiff(type, c("SF", "CH"))
                  }

                  ## calculate distmat if needed
                  if (any(type %in% c("Sil", "D", "COP"))) {
                      if (is.null(a@distmat)) {
                          if (length(a@datalist) == 0L) {
                              warning("Internal CVIs: distmat OR original data needed for indices:",
                                      "\n\tSil\tD\tCOP")

                              type <- setdiff(type, c("Sil", "D", "COP"))

                          } else {
                              distmat <- do.call(a@family@dist,
                                                 args = enlist(x = a@datalist,
                                                               centroids = NULL,
                                                               dots = a@dots))
                          }
                      } else {
                          distmat <- a@distmat
                      }
                  }

                  ## are valid indices still left?
                  if (length(type) == 0L)
                      return(CVIs)

                  ## calculate some values that both Davies-Bouldin indices use
                  if (any(type %in% c("DB", "DBstar"))) {
                      S <- a@clusinfo$av_dist

                      ## distance between centroids
                      distcent <- do.call(a@family@dist,
                                          args = enlist(x = a@centroids,
                                                        centroids = NULL,
                                                        dots = a@dots))
                  }

                  ## calculate global centroids if needed
                  if (any(type %in% c("SF", "CH"))) {
                      N <- length(a@datalist)

                      if (a@type == "partitional") {
                          global_cent <- do.call(a@family@allcent,
                                                 args = enlist(x = a@datalist,
                                                               cl_id = rep(1L, N),
                                                               k = 1L,
                                                               cent = a@datalist[sample(N, 1L)],
                                                               cl_old = rep(0L, N),
                                                               dots = a@dots))
                      } else {
                          global_cent <- a@family@allcent(a@datalist)
                      }

                      dist_global_cent <- do.call(a@family@dist,
                                                  args = enlist(x = a@centroids,
                                                                centroids = global_cent,
                                                                dots = a@dots))

                      dim(dist_global_cent) <- NULL
                  }

                  CVIs <- c(CVIs, sapply(type, function(CVI) {
                      switch(EXPR = CVI,
                             ## Silhouette
                             Sil = {
                                 c_k <- as.numeric(table(a@cluster)[a@cluster])

                                 ab <- lapply(unique(a@cluster), function(k) {
                                     idx <- a@cluster == k

                                     this_a <- rowSums(distmat[idx, idx, drop = FALSE]) / c_k[idx]

                                     this_b <- apply(distmat[idx, !idx, drop = FALSE], 1L, function(row) {
                                         ret <- row / c_k[!idx]
                                         ret <- min(tapply(ret, a@cluster[!idx], sum))
                                         ret
                                     })

                                     data.frame(a = this_a, b = this_b)
                                 })

                                 ab <- do.call(rbind, ab)

                                 sum((ab$b - ab$a) / apply(ab, 1L, max)) / nrow(distmat)
                             },

                             ## Dunn
                             D = {
                                 pairs <- call_pairs(a@k)

                                 deltas <- mapply(pairs[ , 1L], pairs[ , 2L],
                                                  USE.NAMES = FALSE, SIMPLIFY = TRUE,
                                                  FUN = function(i, j) {
                                                      min(distmat[a@cluster == i,
                                                                  a@cluster == j,
                                                                  drop = TRUE])
                                                  })

                                 Deltas <- sapply(1L:a@k, function(k) {
                                     max(distmat[a@cluster == k, a@cluster == k, drop = TRUE])
                                 })

                                 min(deltas) / max(Deltas)
                             },

                             ## Davies-Bouldin
                             DB = {
                                 mean(sapply(1L:a@k, function(k) {
                                     max((S[k] + S[-k]) / distcent[k, -k])
                                 }))
                             },

                             ## Modified DB -> DB*
                             DBstar = {
                                 mean(sapply(1L:a@k, function(k) {
                                     max(S[k] + S[-k]) / min(distcent[k, -k, drop = TRUE])
                                 }))
                             },

                             ## Calinski-Harabasz
                             CH = {
                                 (N - a@k) /
                                     (a@k - 1) *
                                     sum(a@clusinfo$size * dist_global_cent) /
                                     sum(a@cldist[ , 1L, drop = TRUE])
                             },

                             ## Score function
                             SF = {
                                 bcd <- sum(a@clusinfo$size * dist_global_cent) / (N * a@k)
                                 wcd <- sum(a@clusinfo$av_dist)
                                 1 - 1 / exp(exp(bcd - wcd))
                             },

                             ## COP
                             COP = {
                                 1 / nrow(distmat) * sum(sapply(1L:a@k, function(k) {
                                     sum(a@cldist[a@cluster == k, 1L]) / min(apply(distmat[a@cluster != k,
                                                                                           a@cluster == k,
                                                                                           drop = FALSE],
                                                                                   2L,
                                                                                   max))
                                 }))
                             })
                  }))
              }

              CVIs
          })

# ==================================================================================================
# Rand Index from flexclust package
# ==================================================================================================

#' Compare partitions
#'
#' No longer supported directly in \pkg{dtwclust}. Please refer to [cvi()]
#'
#' @name randIndex
#' @rdname randIndex
#' @exportMethod randIndex
#' @keywords internal
#'
#' @param x,y,correct,original See [flexclust::randIndex()].
#'
NULL

#' @rdname randIndex
#' @aliases randIndex,dtwclust,ANY
#' @exportMethod randIndex
#'
setMethod("randIndex", methods::signature(x="dtwclust", y="ANY"),
          function(x, y, correct = TRUE, original = !correct) {
              randIndex(x@cluster, y, correct = correct, original = original)
          })

#' @rdname randIndex
#' @aliases randIndex,ANY,dtwclust
#' @exportMethod randIndex
#'
setMethod("randIndex", methods::signature(x="ANY", y="dtwclust"),
          function(x, y, correct = TRUE, original = !correct) {
              randIndex(x, y@cluster, correct = correct, original = original)
          })

#' @rdname randIndex
#' @aliases randIndex,dtwclust,dtwclust
#' @exportMethod randIndex
#'
setMethod("randIndex", methods::signature(x="dtwclust", y="dtwclust"),
          function(x, y, correct = TRUE, original = !correct) {
              randIndex(x@cluster, y@cluster, correct = correct, original = original)
          })

# ==================================================================================================
# Cluster Similarity from flexclust
# ==================================================================================================

#' Cluster Similarity Matrix
#'
#' No longer supported directly in \pkg{dtwclust}.
#'
#' @name clusterSim
#' @rdname clusterSim
#' @aliases clusterSim,dtwclust-method
#' @exportMethod clusterSim
#' @keywords internal
#'
#' @param object,data,method,symmetric,... See [flexclust::clusterSim()].
#'
#' @seealso
#'
#' [flexclust::clusterSim()]
#'
setMethod("clusterSim", methods::signature(object = "dtwclust"),
          function (object, data = NULL,
                    method = c("shadow", "centers"),
                    symmetric = FALSE, ...)
          {
              method <- match.arg(method)

              if (object@k == 1L)
                  return(matrix(1))

              if (method == "shadow") {
                  if (is.null(data)) {
                      if (!object@control@save.data)
                          stop("No data provided nor found in the dtwclust object.")

                      data <- object@datalist

                  } else {
                      data <- tslist(data)
                  }

                  distmat <- do.call(object@family@dist,
                                     args = enlist(x = data,
                                                   centroids = object@centroids,
                                                   dots = object@dots))

                  r <- t(matrix(apply(distmat, 1L, rank, ties.method = "first"),
                                nrow = ncol(distmat)))
                  cluster <- lapply(1L:2L, function(k) {
                      apply(r, 1L, function(x) which(x == k))
                  })

                  K <- max(cluster[[1L]])
                  z <- matrix(0, ncol = K, nrow = K)
                  for (k in 1L:K) {
                      ok1 <- cluster[[1]] == k
                      if (any(ok1)) {
                          for (n in 1L:K) {
                              if (k != n) {
                                  ok2 <- ok1 && cluster[[2L]] == n
                                  if (any(ok2)) {
                                      z[k, n] <- 2 * sum(distmat[ok2, k] / (distmat[ok2, k] + distmat[ok2, n]))
                                  }
                              }
                          }

                          z[k, ] <- z[k, ] / sum(ok1)
                      }
                  }

                  diag(z) <- 1

                  if (symmetric)
                      z <- (z + t(z)) / 2

              } else {
                  z <- do.call(object@family@dist,
                               args = enlist(x = object@centroids,
                                             centroids = object@centroids,
                                             dots = object@dots))

                  z <- 1 - z / max(z)
              }

              z
          })

# ==================================================================================================
# Validity and coercion methods for control
# ==================================================================================================

setValidity("dtwclustControl",
            function(object) {
                if (!is.null(object@window.size) && object@window.size < 1L)
                    return("Window size must be positive if provided")

                if (!(object@norm %in% c("L1", "L2")))
                    return("norm must be either L1 or L2")

                if (object@dba.iter < 0L)
                    return("DBA iterations must be positive")

                if (object@iter.max < 0L)
                    return("Maximum iterations must be positive")

                if (object@nrep < 1L)
                    return("Number of repetitions must be at least one")

                if (object@fuzziness <= 1)
                    return("Fuzziness exponent should be greater than one")

                if (object@delta < 0)
                    return("Delta should be positive")

                for (sl in slotNames("dtwclustControl")) {
                    if (sl == "packages") next

                    if (length(slot(object, sl)) > 1L)
                        return(paste0("Control parameters should only have one element ",
                                      "(", sl, " had too many)"))
                }

                TRUE
            })

setAs("list", "dtwclustControl",
      function(from, to) {
          ctrl <- methods::new(to)

          num <- c("delta", "fuzziness")

          for (arg in names(from)) {
              val <- from[[arg]]

              if (is.numeric(val) && !(arg %in% num))
                  val <- as.integer(val)

              methods::slot(ctrl, arg) <- val
          }

          methods::validObject(ctrl)

          ctrl
      })

setAs("NULL", "dtwclustControl",
      function(from, to) {
          methods::new(to)
      })

# ==================================================================================================
# Functions to support package 'clue'
# ==================================================================================================

#' @method n_of_classes dtwclust
#' @export
#'
n_of_classes.dtwclust <- function(x) {
    x@k
}

#' @method n_of_objects dtwclust
#' @export
#'
n_of_objects.dtwclust <- function(x) {
    length(x@cluster)
}

#' @method cl_class_ids dtwclust
#' @export
#'
cl_class_ids.dtwclust <- function(x) {
    clue::as.cl_class_ids(x@cluster)
}

#' @method as.cl_membership dtwclust
#' @export
#'
as.cl_membership.dtwclust <- function(x) {
    as.cl_membership(x@cluster)
}

#' @method cl_membership dtwclust
#' @export
#'
cl_membership.dtwclust <- function(x, k = n_of_classes(x)) {
    as.cl_membership(x)
}

#' @method is.cl_partition dtwclust
#' @export
#'
is.cl_partition.dtwclust <- function(x) {
    TRUE
}

#' @method is.cl_hard_partition dtwclust
#' @export
#'
is.cl_hard_partition.dtwclust <- function(x) {
    x@type != "fuzzy"
}

#' @method is.cl_hierarchy dtwclust
#' @export
#'
is.cl_hierarchy.dtwclust <- function(x) {
    x@type == "hierarchical"
}

#' @method is.cl_dendrogram dtwclust
#' @export
#'
is.cl_dendrogram.dtwclust <- function(x) {
    x@type == "hierarchical"
}
