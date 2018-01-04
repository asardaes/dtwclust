# ==================================================================================================
# Internal global variables
# ==================================================================================================

supported_clusterings <- c("partitional", "hierarchical", "fuzzy", "tadpole")

distances_known <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd", "dtw_basic", "gak", "sdtw")
distances_included <- c("dtw_lb", "lb_keogh", "lb_improved", "sbd", "dtw_basic", "gak", "sdtw")
distances_difflength <- c("dtw", "dtw2", "sbd", "dtw_basic", "gak", "sdtw")
distances_multivariate <- c("dtw", "dtw2", "dtw_basic", "gak", "sdtw")

centroids_included <- c("mean", "median", "shape", "dba", "pam", "fcm", "fcmdd", "sdtw_cent")
centroids_fuzzy <- c("fcm", "fcmdd")
centroids_nonfuzzy <- setdiff(centroids_included, centroids_fuzzy)
centroids_difflength <- c("dba", "pam", "shape", "fcmdd", "sdtw_cent")

control_classes <- c(partitional = "PtCtrl",
                     hierarchical = "HcCtrl",
                     fuzzy = "FzCtrl",
                     tadpole = "TpCtrl",
                     args = "TscArgs")
