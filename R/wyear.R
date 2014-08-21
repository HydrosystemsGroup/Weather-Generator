#' Wavelet analysis of annual timeseries
#'
#' @param x an array of dates in formats accepted by lubridate
#' @param start.month initial month of the water year (default is 10 for October)
#' @export
#' @examples
#' library(lurbidate)
#' wyear(ymd(c("2000-01-01", "2000-09-30", "2000-10-01", "2000-12-31")), start.month=10)
wyear <- function(x, start.month=10) {
  x.wyr <- ifelse(lubridate::month(x)>=start.month, lubridate::year(x)+1, lubridate::year(x))
  return(x.wyr)
}