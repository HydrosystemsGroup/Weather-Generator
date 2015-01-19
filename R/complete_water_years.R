#' Get Indices of Complete Water Years
#'
#' @param x vector of dates
#' @param start_month first month of the water year
#' @export
#' @examples
#' dates <- seq.Date(from=as.Date("2000-01-01"), to=as.Date("2004-12-31"), by="1 day")
#' dates_complete <- dates[complete_water_years(dates, start_month=10)] # 2000-10-01 -> 20014-09-30
complete_water_years <- function(x, start_month=10) {
  wyears <- wyear(x, start_month = start_month)
  wyears_counts <- aggregate(wyears, by=list(WYEAR=wyears), FUN=length)
  wyears_complete <- wyears_counts[which(wyears_counts[['x']] >= 365), 'WYEAR']
  which(wyears %in% wyears_complete)
}
