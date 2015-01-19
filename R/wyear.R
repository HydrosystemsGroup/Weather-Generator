#' Get water year component of a date-time
#'
#' @param x a vector of date-times
#' @param start_month initial month of the water year (default is 10 for October)
#' @export
#' @examples
#' wyear(lubridate::ymd("2000-10-01"), start_month=10)
#' wyear(lubridate::ymd(c("2000-01-01", "2000-09-30", "2000-10-01", "2000-12-31")), start_month=10)
wyear <- function(x, start_month=10) {
  x.wyr <- ifelse(lubridate::month(x)>=start_month, lubridate::year(x)+1, lubridate::year(x))
  return(x.wyr)
}
