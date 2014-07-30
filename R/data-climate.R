#' Climate dataset
#'
#' This dataset contains daily climate data provided by Scott Steinschneider.
#' It appears to be GHCND-Daily data for a number of stations, with gaps removed.
#' Some files end in "MANUAL" so there may have been some manual adjustments.
#' It will do for now for testing.
#'
#' @section Dimensions:
#'
#' \itemize{
#'   \item \code{FILE}: name of original file
#'   \item \code{DATE}: dates as Date objects
#' }
#'
#' @section Measures:
#'
#' \itemize{
#'   \item \code{PRCP}: daily precipitation (mm/day)
#'   \item \code{TMIN}, \code{TMAX}: daily minimum and maximum temperature (degC)
#' }
#' @docType data
#' @name climate
#' @usage climate
#' @format A data.frame with 1,027,540 observations.
#' @examples
#' climate
NULL