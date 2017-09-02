#' Subset of character trajectories data set
#'
#' Subset: only 5 examples of each considered character. See details.
#'
#' @name uciCT
#' @aliases ucict CharTraj CharTrajLabels CharTrajMV
#'
#' @format
#'
#' Lists with 100 elements. Each element is a time series. Labels included as factor vector.
#'
#' @details
#'
#' Quoting the source:
#'
#' "Multiple, labelled samples of pen tip trajectories recorded whilst writing individual
#' characters. All samples are from the same writer, for the purposes of primitive extraction. Only
#' characters with a single pen-down segment were considered."
#'
#' The subset included in `CharTraj` has only 5 examples of the X velocity for each character. A
#' vector with labels is also loaded in `CharTrajLabels`.
#'
#' The subset included in `CharTrajMV` has 5 examples too, but includes tip force as well as X and Y
#' velocity. Each element of the list is a multivariate series with 3 variables.
#'
#' Please note that even though both `CharTraj` and `CharTrajMV` have the same series names, the
#' actual series in each subset are **not** the same, i.e., `CharTraj$A.V1` is not in
#' `CharTrajMV$A.V1`.
#'
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Character+Trajectories}
#'
NULL

#' Results of timing experiments
#'
#' This is the list with data frames containing the results of the timing experiments vignette
#' included with \pkg{dtwclust}. See `browseVignettes("dtwclust")`.
#'
#' @name dtwclustTimings
#'
#' @format
#'
#' The results are organized into different data frames and saved in one list with nested lists.
#' For more details, refer to the included vignette or the scripts available at
#' \url{https://github.com/asardaes/dtwclust/tree/master/timing-experiments}.
#'
#' @source Refer to the timing experiments vignette.
#'
NULL
