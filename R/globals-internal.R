# ==================================================================================================
# Internal global variables
# ==================================================================================================

distances_included <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd", "dtw_basic", "gak")
distances_difflength <- c("dtw", "dtw2", "sbd", "dtw_basic", "gak")
distances_multivariate <- c("dtw", "dtw2", "dtw_basic", "gak")

centroids_included <- c("mean", "median", "shape", "dba", "pam", "fcm", "fcmdd")
centroids_fuzzy <- c("fcm", "fcmdd")
centroids_difflength <- c("dba", "pam", "shape", "fcmdd")
