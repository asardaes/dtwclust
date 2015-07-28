#' Subset of character trajectories data set
#'
#' Subset: only 5 examples of X velocity. See details.
#'
#' Quoting the source:
#'
#' "Multiple, labelled samples of pen tip trajectories recorded whilst writing individual characters.
#' All samples are from the same writer, for the purposes of primitive extraction. Only characters with a
#' single pen-down segment were considered."
#'
#' The subset included here (\code{CharTraj}) has only 5 examples of the X velocity for each character.
#' A vector with labels is also loaded in \code{CharTrajLabels}.
#'
#' @name uciCT
#'
#' @format A list with 100 elements. Each element is a time series. Labels included as factor vector.
#'
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Character+Trajectories}
#'
#' @aliases ucict CharTraj CharTrajLabels
#'
NULL
