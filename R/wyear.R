#' Get water year component of a date-time
#'
#' @param x a vector of date-times
#' @param start.month initial month of the water year (default is 10 for October)
#' @export
#' @examples
#' wyear(lubridate::ymd("2000-10-01"), start.month=10)
#' wyear(lubridate::ymd(c("2000-01-01", "2000-09-30", "2000-10-01", "2000-12-31")), start.month=10)
wyear <- function(x, start.month=10) {
  x.wyr <- ifelse(lubridate::month(x)>=start.month, lubridate::year(x)+1, lubridate::year(x))
  return(x.wyr)
}

#' Get julian day component by water year of a date-time (aka water day)
#'
#' @param x a vector of date-times
#' @param start.month initial month of the water year (default is 10 for October)
#' @export
#' @examples
#' waterday(lubridate::ymd("2000-10-01"), start.month=10)
#' waterday(lubridate::ymd(c("2000-01-01", "2000-09-30", "2000-10-01", "2000-12-31")), start.month=10)
waterday <- function(x, start.month=10) {
  start_jday <- lubridate::yday(lubridate::ymd(paste(2001, start.month, 1, sep='-'))) - 1
  wday <- lubridate::yday(x)-(start_jday+lubridate::leap_year(x))
  wday <- ifelse(wday<=0, 365+wday+leap_year(x), wday)
  wday
}

#' Get a range of water days around some center value
#'
#' @param x central water day
#' @param n length of range
#' @export
#' @examples
#' waterday_range(x=1, n=7)
#' waterday_range(x=10, n=7)
#' waterday_range(x=363, n=7)
waterday_range <- function(x, n=7) {
  if (n %% 2 == 1) {
    l <- (n-1)/2
    u <- (n-1)/2
  } else {
    l <- n/2
    u <- n/2 - 1
  }
  rng <- seq(x-l, x+u)
  rng <- ifelse(rng<=0, 365+rng, rng)
  rng <- ifelse(rng>365, rng-365, rng)
  rng
}
