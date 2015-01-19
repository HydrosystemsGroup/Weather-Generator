#' Create Sequence of Dates by Water Year
#'
#' @param n_year number of simulation years
#' @param start_month first month of the water year
#' @param start_water_year initial water year (e.g. if start_water_year=2000 and start_month=10 then the first date is 1999-10-01)
#' @param include_leap_days boolean flag to include leap days (default=FALSE)
#' @export
#' @examples
#' wyear_date_range(n_year=10, start_month=10, start_water_year=2000, include_leap_days=FALSE)
wyear_date_range <- function(n_year, start_month=10, start_water_year=2000, include_leap_days=FALSE) {
  if (include_leap_days) {
    stop('Leap days not yet supported')
  }

  if (start_month != 1) {
    start_year <- start_water_year - 1
  } else {
    start_year <- start_water_year
  }

  dates <- seq.Date(from=as.Date(paste(start_year, start_month, 1, sep='-')),
                    to=as.Date(paste(start_year+n_year, start_month, 1, sep='-')) - lubridate::days(1),
                    by=1)

  # remove leap days
  leap_days <- which(lubridate::month(dates)==2 & lubridate::day(dates)==29)
  if (length(leap_days) > 0) {
    dates <- dates[-leap_days]
  }

  stopifnot(length(dates) %% 365 == 0)

  dates
}
