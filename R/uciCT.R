#' Subset of character trajectories data set
#'
#' Subset: only 5 examples of each considered character. See details.
#'
#' @name uciCT
#' @aliases ucict CharTraj CharTrajLabels CharTrajMV
#'
#' @format
#'
#' Lists with 100 elements each. Each element is a time series. Labels included as factor vector.
#'
#' @details
#'
#' Quoting the source:
#'
#' "Multiple, labelled samples of pen tip trajectories recorded whilst writing individual
#' characters. All samples are from the same writer, for the purposes of primitive extraction. Only
#' characters with a single pen-down segment were considered."
#'
#' The subset included in \code{CharTraj} has only 5 examples of the X velocity for each character.
#' A vector with labels is also loaded in \code{CharTrajLabels}.
#'
#' The subset included in \code{CharTrajMV} has 5 examples too, but includes tip force as well as X
#' and Y velocity. Each element of the list is a multivariate series with 3 variables.
#'
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Character+Trajectories}
#'
NULL
