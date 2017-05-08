# ========================================================================================================
# Custom functions to calculate centroids
# ========================================================================================================

all_cent <- function(case = NULL, distmat = NULL, distfun, control, fuzzy = FALSE) {
    ## --------------------------------------------------------------------------------------------
    ## pam
    pam_cent <- function(x, x_split, cl_id, id_changed, ...) {
        if (is.null(distmat)) {
            new_cent <- lapply(x_split, function(xsub) {
                distmat <- distfun(xsub, xsub)

                d <- apply(distmat, 1L, sum)

                i_cent <- xsub[[which.min(d)]]

                ## update indices, aggregated at the end of main allcent function
                attr(i_cent, "id_cent") <- pmatch(names(xsub[which.min(d)]), names(x))

                i_cent
            })

        } else {
            id_x <- lapply(id_changed, function(cl_num) which(cl_id == cl_num))

            new_cent <- lapply(id_x, function(i_x) {
                d <- apply(distmat[i_x, i_x, drop = FALSE], 1L, sum)

                i_cent <- x[[i_x[which.min(d)]]]

                ## update indices, aggregated at the end of main allcent function
                attr(i_cent, "id_cent") <- i_x[which.min(d)]

                i_cent
            })
        }

        ## return
        new_cent
    }

    ## --------------------------------------------------------------------------------------------
    ## shape
    shape_cent <- function(x_split, cent, ...) {
        x_split <- split_parallel(x_split)
        cent <- split_parallel(cent)

        new_cent <- foreach(x_split = x_split,
                            cent = cent,
                            .combine = c,
                            .multicombine = TRUE,
                            .packages = "dtwclust") %op% {
                                mapply(x_split, cent,
                                       SIMPLIFY = FALSE,
                                       FUN = function(x, c) {
                                           new_c <- shape_extraction(x, c)

                                           ## return
                                           new_c
                                       })
                            }

        ## return
        new_cent
    }

    ## --------------------------------------------------------------------------------------------
    ## dba
    dba_cent <- function(x, x_split, cent, id_changed, cl_id, ...) {
        ## not all arguments are used, but I want them to be isolated from ...
        dots <- list(...)

        x_split <- split_parallel(x_split)
        cent <- split_parallel(cent)

        new_cent <- foreach(x_split = x_split,
                            cent = cent,
                            .combine = c,
                            .multicombine = TRUE,
                            .packages = "dtwclust",
                            .export = c("control", "enlist")) %op% {
                                mapply(x_split, cent,
                                       SIMPLIFY = FALSE,
                                       FUN = function(x, c) {
                                           new_c <- do.call(DBA,
                                                            enlist(X = x,
                                                                   centroid = c,
                                                                   norm = control@norm,
                                                                   window.size = control@window.size,
                                                                   max.iter = control@dba.iter,
                                                                   delta = control@delta,
                                                                   error.check = FALSE,
                                                                   dots = dots))

                                           ## return
                                           new_c
                                       })
                            }

        ## return
        new_cent
    }

    ## --------------------------------------------------------------------------------------------
    ## mean
    mean_cent <- function(x_split, ...) {
        new_cent <- lapply(x_split, function(xx) {
            if (is_multivariate(xx)) {
                ## multivariate
                ncols <- ncol(xx[[1L]]) # number of dimensions should be equal
                ncols <- rep(1L:ncols, length(xx))

                xx <- do.call(cbind, xx)
                xx <- split.data.frame(t(xx), ncols)
                do.call(cbind, lapply(xx, colMeans))

            } else {
                ## univariate
                xx <- do.call(rbind, xx)
                colMeans(xx)
            }
        })

        ## return
        new_cent
    }

    ## --------------------------------------------------------------------------------------------
    ## median
    median_cent <- function(x_split, ...) {
        new_cent <- lapply(x_split, function(xx) {
            if (is_multivariate(xx)) {
                ## multivariate
                ncols <- ncol(xx[[1L]]) # number of dimensions should be equal
                ncols <- rep(1L:ncols, length(xx))

                xx <- do.call(cbind, xx)
                xx <- split.data.frame(t(xx), ncols)
                do.call(cbind, lapply(xx, colMedians))

            } else {
                ## univariate
                xx <- do.call(rbind, xx)
                colMedians(xx)
            }
        })

        ## return
        new_cent
    }

    ## --------------------------------------------------------------------------------------------
    ## fcm
    fcm_cent <- function(x, u, k, ...) {
        ## utils.R
        if (is_multivariate(x)) {
            mv <- reshape_multivariate(x, NULL)

            cent <- lapply(mv$series, fcm_cent, u = u)
            cent <- lapply(1L:k, function(idc) {
                sapply(cent, function(c) { c[idc, , drop = TRUE] })
            })
            cent <- lapply(cent, "dimnames<-", NULL)

            cent

        } else {
            cent <- t(u) %*% do.call(rbind, x)
            apply(cent, 2L, "/", e2 = colSums(u))
        }
    }

    ## --------------------------------------------------------------------------------------------
    ## fcmdd
    fcmdd_cent <- function(x, u, k, ...) {
        q <- distmat %*% u
        idc <- apply(q, 2L, which.min)

        cent <- x[idc]
        attr(cent, "id_cent") <- idc

        cent
    }

    ## --------------------------------------------------------------------------------------------
    ## allcent
    if (fuzzy) {
        ## function created here to capture objects of this environment (closure)
        allcent <- function(x, cl_id, k, cent, cl_old, ...) {
            ## cent and cl_old are unused here, but R complains if signatures don't match
            u <- cl_id ^ control@fuzziness

            cent <- do.call(paste0(case, "_cent"),
                            enlist(x = x,
                                   u = u,
                                   k = k,
                                   dots = list(...)))

            ## coerce back to list
            any2list(cent)
        }
    } else {
        allcent <- function(x, cl_id, k, cent, cl_old, ...) {
            ## Check which clusters changed
            if (all(cl_old == 0L)) {
                id_changed <- sort(unique(cl_id))

            } else {
                id_changed <- cl_id != cl_old
                id_changed <- union(cl_id[id_changed], cl_old[id_changed])
            }

            if (length(id_changed) == 0L) {
                return(cent)
            }

            ## Split data according to cluster memebership
            x_split <- split(x, factor(cl_id, levels = 1L:k))

            ## In case of empty new clusters
            empty_clusters <- which(lengths(x_split) == 0L)
            id_changed <- setdiff(id_changed, empty_clusters)

            ## Calculate new centroids
            new_cent <- do.call(paste0(case, "_cent"),
                                enlist(x = x,
                                       x_split = x_split[id_changed],
                                       cent = cent[id_changed],
                                       id_changed = id_changed,
                                       cl_id = cl_id,
                                       dots = list(...)))

            cent[id_changed] <- new_cent

            ## Any empty clusters?
            num_empty <- length(empty_clusters)

            ## If so, initialize new clusters
            if (num_empty > 0L) {
                ## Make sure no centroid is repeated (especially in case of PAM)
                any_rep <- logical(num_empty)

                while(TRUE) {
                    id_cent_extra <- sample(length(x), num_empty)
                    extra_cent <- x[id_cent_extra]

                    for (id_extra in 1L:num_empty) {
                        any_rep[id_extra] <- any(sapply(cent, function(i.centroid) {
                            if (length(i.centroid) != length(extra_cent[[id_extra]]))
                                ret <- FALSE
                            else
                                ret <- all(i.centroid == extra_cent[[id_extra]])

                            ret
                        }))

                        if (case == "pam")
                            attr(extra_cent[[id_extra]], "id_cent") <- id_cent_extra[id_extra]
                    }

                    if (all(!any_rep))
                        break
                }

                cent[empty_clusters] <- extra_cent
            }

            ## aggregate updated indices
            if (case == "pam")
                attr(cent, "id_cent") <- sapply(cent, attr, which = "id_cent")

            cent
        }
    }

    allcent
}
