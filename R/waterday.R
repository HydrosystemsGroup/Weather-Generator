#' Get julian day component by water year of a date-time (aka water day)
#'
#' @param x a vector of date-times
#' @param start_month initial month of the water year (default is 10 for October)
#' @export
#' @examples
#' waterday(lubridate::ymd("2000-10-01"), start_month=10)
#' waterday(lubridate::ymd(c("2000-01-01", "2000-09-30", "2000-10-01", "2000-12-31")), start_month=10)
waterday <- function(x, start_month=10) {
  start_jday <- lubridate::yday(lubridate::ymd(paste(2001, start_month, 1, sep='-'))) - 1
  wday <- lubridate::yday(x)-(start_jday+lubridate::leap_year(x))
  wday <- ifelse(wday<=0, 365+wday+lubridate::leap_year(x), wday)
  wday
}