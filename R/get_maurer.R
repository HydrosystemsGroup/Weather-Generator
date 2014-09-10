#' Retrieve monthly timeseries for given location from maurer dataset
#'
#' @param lat latitude in decimal degrees North
#' @param lon longitude in decimal degrees East (negative for West)
#' @export
#' @examples
#' get_maurer_mon(42, -71)
#'
get_maurer_mon <- function(lat, lon) {
  data(maurer)
  distances <- sqrt((maurer$grid$LAT-lat)^2 + (maurer$grid$LON-lon)^2)
  pos <- c(maurer$grid[which.min(distances),])
  maurer_ts <- subset(maurer$data, LAT==pos$LAT & LON==pos$LON)
  return(maurer_ts)
}

#' Aggregate monthly timeseries from maurer to annual timesteps
#'
#' @param x data frame containing monthly maurer timeseries
#' @export
#' @examples
#' aggregate_maurer_mon(get_maurer_mon(42, -71))
#'
aggregate_maurer_mon <- function(x) {
  require(dplyr)
  x <- dplyr::mutate(x, N_DAY=lubridate::days_in_month(DATE))
  grp <- dplyr::group_by(x, LAT, LON, YEAR)
  x.yr <- dplyr::summarise(grp,
                           PRCP=sum(PRCP),
                           TMAX=sum(TMAX*N_DAY)/sum(N_DAY),
                           TMIN=sum(TMIN*N_DAY)/sum(N_DAY),
                           WIND=sum(WIND*N_DAY)/sum(N_DAY))
  x.yr <- dplyr::select(x.yr, YEAR, PRCP, TMAX, TMIN, WIND, LAT, LON)
  return(x.yr)
}