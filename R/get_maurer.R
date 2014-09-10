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
  return(list(location=list(latitude=lat, longitude=lon),
              nearest=list(latitude=pos$LAT, longitude=pos$LON),
              data=maurer_ts))
}