# ==================================================================================================
# Internal global variables
# ==================================================================================================

supported_clusterings <- c("partitional", "hierarchical", "fuzzy", "tadpole")

distances_known <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd", "dtw_basic", "gak")
distances_included <- c("dtw_lb", "lbk", "lbi", "sbd", "dtw_basic", "gak")
distances_difflength <- c("dtw", "dtw2", "sbd", "dtw_basic", "gak")
distances_multivariate <- c("dtw", "dtw2", "dtw_basic", "gak")

centroids_included <- c("mean", "median", "shape", "dba", "pam", "fcm", "fcmdd")
centroids_fuzzy <- c("fcm", "fcmdd")
centroids_nonfuzzy <- setdiff(centroids_included, centroids_fuzzy)
centroids_difflength <- c("dba", "pam", "shape", "fcmdd")

control_classes <- c(partitional = "PtCtrl",
                     hierarchical = "HcCtrl",
                     fuzzy = "FzCtrl",
                     tadpole = "TpCtrl",
                     args = "TscArgs")
