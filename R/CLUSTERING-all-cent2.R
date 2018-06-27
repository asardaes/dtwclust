# ==================================================================================================
# Helpers
# ==================================================================================================

# column-wise medians
#' @importFrom stats median
#'
colMedians <- function(mat) { apply(mat, 2L, stats::median) }

# Reinitialize empty clusters
reinit_clusters <- function(x, cent, cent_case, num_empty, empty_clusters, control) {
    # Make sure no centroid is repeated (especially in case of PAM)
    any_rep <- logical(num_empty)
    while(TRUE) {
        id_cent_extra <- sample(length(x), num_empty)
        extra_cent <- x[id_cent_extra]
        for (id_extra in 1L:num_empty) {
            any_rep[id_extra] <- any(sapply(cent, function(i.centroid) {
                identical(i.centroid, extra_cent[[id_extra]])
            }))
            if (cent_case == "pam")
                control$distmat$id_cent[empty_clusters[id_extra]] <- id_cent_extra[id_extra]
        }
        if (all(!any_rep)) break
    }
    extra_cent
}

# ==================================================================================================
# Custom functions to calculate centroids
# ==================================================================================================

#' @importFrom Matrix rowSums
#'
all_cent2 <- function(case = NULL, control) {
    # ----------------------------------------------------------------------------------------------
    # pam
    pam_cent <- function(x, x_split, cent, id_changed, cl_id, ...) {
        id_x <- lapply(id_changed, function(cl_num) { which(cl_id == cl_num) })
        # return
        Map(id_x, id_changed, f = function(i_x, i_cl) {
            # must specify because the one defined in this environment has same name (historic)
            dtwclust::pam_cent(x, ids = i_x, distmat = control$distmat, ...,
                               .i_cl_ = i_cl, error.check = FALSE)
        })
    }

    # ----------------------------------------------------------------------------------------------
    # shape
    shape_cent <- function(x, x_split, cent, id_changed, cl_id, ...) {
        # not all arguments are used, but I want them to be isolated from ellipsis
        dots <- list(...)
        dots$error.check <- FALSE
        x_split <- split_parallel(x_split)
        cent <- split_parallel(cent)
        # return
        foreach(x_split = x_split, cent = cent,
                .combine = c,
                .multicombine = TRUE,
                .packages = "dtwclust",
                .export = c("enlist")) %op% {
                    Map(x_split, cent, f = function(x, cent) {
                        do.call(shape_extraction, enlist(X = x, centroid = cent, dots = dots), TRUE)
                    })
                }
    }

    # ----------------------------------------------------------------------------------------------
    # dba
    dba_cent <- function(x, x_split, cent, id_changed, cl_id, ...) {
        # not all arguments are used, but I want them to be isolated from ellipsis
        dots <- list(...)
        dots$error.check <- FALSE
        x_split <- split_parallel(x_split)
        cent <- split_parallel(cent)
        # return
        foreach(x_split = x_split, cent = cent,
                .combine = c,
                .multicombine = TRUE,
                .packages = "dtwclust",
                .export = c("enlist")) %op% {
                    Map(x_split, cent, f = function(x, cent) {
                        do.call(DBA, enlist(X = x, centroid = cent, dots = dots), TRUE)
                    })
                }
    }

    # ----------------------------------------------------------------------------------------------
    # sdtw_cent
    sdtw_cent_cent <- function(x, x_split, cent, id_changed, cl_id, ...) {
        # not all arguments are used, but I want them to be isolated from ellipsis
        dots <- list(...)
        dots$error.check <- FALSE
        x_split <- split_parallel(x_split)
        cent <- split_parallel(cent)
        # return
        foreach(x_split = x_split, cent = cent,
                .combine = c,
                .multicombine = TRUE,
                .packages = "dtwclust",
                .export = c("enlist")) %op% {
                    Map(x_split, cent, f = function(x, cent) {
                        do.call(sdtw_cent, enlist(series = x, centroid = cent, dots = dots), TRUE)
                    })
                }
    }

    # ----------------------------------------------------------------------------------------------
    # mean
    mean_cent <- function(x_split, ...) {
        # return
        lapply(x_split, function(xx) {
            if (is_multivariate(xx)) {
                ncols <- ncol(xx[[1L]]) # number of dimensions should be equal
                ncols <- rep(1L:ncols, length(xx))
                xx <- do.call(cbind, xx, TRUE)
                xx <- split.data.frame(t(xx), ncols)
                do.call(cbind, lapply(xx, colMeans), TRUE)
            }
            else {
                xx <- do.call(cbind, xx, TRUE)
                rowMeans(xx)
            }
        })
    }

    # ----------------------------------------------------------------------------------------------
    # median
    median_cent <- function(x_split, ...) {
        # return
        lapply(x_split, function(xx) {
            if (is_multivariate(xx)) {
                ncols <- ncol(xx[[1L]]) # number of dimensions should be equal
                ncols <- rep(1L:ncols, length(xx))
                xx <- do.call(cbind, xx, TRUE)
                xx <- split.data.frame(t(xx), ncols)
                do.call(cbind, lapply(xx, colMedians), TRUE)
            }
            else {
                xx <- do.call(rbind, xx, TRUE)
                colMedians(xx)
            }
        })
    }

    # ----------------------------------------------------------------------------------------------
    # fcm
    fcm_cent <- function(x, u, k, ...) {
        # UTILS-utils.R
        if (is_multivariate(x)) {
            mv <- reshape_multivariate(x, NULL)
            cent <- lapply(mv$series, fcm_cent, u = u)
            cent <- lapply(1L:k, function(idc) {
                sapply(cent, function(cent) { cent[idc, , drop = TRUE] })
            })
            cent <- lapply(cent, "dimnames<-", NULL)
            cent
        }
        else {
            cent <- t(u) %*% do.call(rbind, x, TRUE)
            apply(cent, 2L, "/", e2 = colSums(u))
        }
    }

    # ----------------------------------------------------------------------------------------------
    # fcmdd
    fcmdd_cent <- function(x, u, k, ...) {
        q <- control$distmat$distmat %*% u
        idc <- apply(q, 2L, which.min)
        cent <- x[idc]
        control$distmat$id_cent <- idc
        cent
    }

    # ----------------------------------------------------------------------------------------------
    # allcent
    if (case %in% c("fcm", "fcmdd")) {
        # function created here to capture objects of this environment (closure)
        allcent <- function(x, cl_id, k, cent, cl_old, ...) {
            # cent and cl_old are unused here, but R complains if signatures don't match
            x <- tslist(x)
            u <- cl_id ^ control$fuzziness
            cent <- do.call(paste0(case, "_cent"),
                            enlist(x = x,
                                   u = u,
                                   k = k,
                                   dots = list(...)),
                            TRUE)
            # coerce back to list
            tslist(cent)
        }
    }
    else {
        allcent <- function(x, cl_id, k, cent, cl_old, ...) {
            x <- tslist(x)
            cent <- tslist(cent)
            # Check which clusters changed
            if (all(cl_old == 0L)) {
                id_changed <- sort(unique(cl_id))
            }
            else {
                id_changed <- cl_id != cl_old
                id_changed <- union(cl_id[id_changed], cl_old[id_changed])
            }
            if (length(id_changed) == 0L) return(cent)

            # Split data according to cluster memebership
            x_split <- split(x, factor(cl_id, levels = 1L:k))
            # In case of empty new clusters
            empty_clusters <- which(lengths(x_split) == 0L)
            id_changed <- setdiff(id_changed, empty_clusters)
            # Calculate new centroids
            new_cent <- do.call(paste0(case, "_cent"),
                                enlist(x = x,
                                       x_split = x_split[id_changed],
                                       cent = cent[id_changed],
                                       id_changed = id_changed,
                                       cl_id = cl_id,
                                       dots = list(...)),
                                TRUE)
            cent[id_changed] <- new_cent

            # Any empty clusters?
            num_empty <- length(empty_clusters)
            # If so, initialize new clusters
            if (num_empty > 0L)
                cent[empty_clusters] <- reinit_clusters(x, cent, case, num_empty,
                                                        empty_clusters, control)
            # return
            cent
        }
    }
    # return closure
    allcent
}
