#' Climate cities dataset
#'
#' This dataset contains daily climate data for three cities in the northeast
#' from the Maurer dataset. The cities include Boston, Providence, and Worcester.
#' The dataset is provided as list of data frames. The dataset was generated
#' using the R script: /scripts/postgres/climate_cities.R
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
#'   \item \code{WIND}: daily mean windspeed (m/s)
#' }
#' @docType data
#' @name climate_cities
#' @usage climate_cities
#' @format A list of length three, each element a data frame with 22,645 rows
#' @examples
#' climate_cities[['boston']]
NULL