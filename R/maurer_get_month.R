#' Retrieve monthly timeseries for given location from maurer dataset
#'
#' @param lat latitude in decimal degrees North
#' @param lon longitude in decimal degrees East (negative for West)
#' @export
#' @examples
#' maurer_get_month(42, -71)
#'
maurer_get_month <- function(lat, lon) {
  data(maurer)
  distances <- sqrt((maurer$grid$LAT-lat)^2 + (maurer$grid$LON-lon)^2)
  pos <- c(maurer$grid[which.min(distances),])
  maurer_ts <- subset(maurer$data, LAT==pos$LAT & LON==pos$LON)
  maurer_ts <- dplyr::mutate(maurer_ts, TAVG=(TMIN+TMAX)/2)
  return(maurer_ts)
}
