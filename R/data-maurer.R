#' Maurer meterological dataset
#'
#' This dataset contains a subset of the monthly meteorological datasets created
#' by Maurer et al. (2002). Only locations in the northeast region with latitudes
#' above 40 N are included. This dataset was automatically generated
#' using the makefile and r scripts located in the \code{/scripts/maurer/}
#' directory.
#'
#' The object \code{maurer} is a two-element list containing
#'
#' \itemize{
#'   \item \code{data}: data frame of monthly climate timeseries
#'   \item \code{grid}: data frame containing the file name, latitude and longitude of the locations used in the regional average
#' }
#'
#' @section Measures:
#'
#' The \code{data} data frame contains the following variables.
#'
#' \itemize{
#'   \item \code{PRCP}: precipitation (mm)
#'   \item \code{TMIN}, \code{TMAX}: minimum and maximum temperature (degC)
#'   \item \code{WIND}: wind speed (m/s)
#' }
#' @docType data
#' @name maurer
#' @usage maurer
#' @references
#'   Maurer, E.P., A.W. Wood, J.C. Adam, D.P. Lettenmaier, and B. Nijssen, 2002,
#'   A Long-Term Hydrologically-Based Data Set of Land Surface Fluxes and States
#'   for the Conterminous United States, J. Climate 15, 3237-3251.
#'   \url{http://www.engr.scu.edu/~emaurer/gridded_obs/index_gridded_obs.html}
#' @examples
#' maurer
NULL