# ========================================================================================================
# Membership update for fuzzy c-means clustering
# ========================================================================================================

fcm_cluster <- function(distmat, m) {
    cprime <- apply(distmat, 1L, function(dist_row) { sum( (1 / dist_row) ^ (2 / (m - 1)) ) })
    u <- 1 / apply(distmat, 2L, function(dist_col) { cprime * dist_col ^ (2 / (m - 1)) })
    if (is.null(dim(u))) u <- rbind(u) # for predict generic
    u[is.nan(u)] <- 1 # in case fcmdd is used
    u
}

# ========================================================================================================
# Fuzzy objective function
# ========================================================================================================

fuzzy_objective <- function(u, distmat, m) { sum(u^m * distmat^2) }
