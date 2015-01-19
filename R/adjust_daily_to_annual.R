#' Adjust Daily Values to Match Annual Sums
#'
#' @param x_day vector of daily values
#' @param years_day vector of years corresponding to x_day
#' @param x_yr vector of annual sums
#' @param years_yr vector of years corresponding to x_yr
#' @param min_ratio minimum adjustment ratio
#' @param max_ratio maximum adjustment ratio
#' @export
adjust_daily_to_annual <- function(x_day, years_day, x_yr, years_yr, min_ratio=0.9, max_ratio=1.1) {
  stopifnot(length(x_day)==length(years_day))
  stopifnot(length(x_yr)==length(years_yr))

  # create data frames
  daily <- data.frame(YEAR=years_day,
                      VALUE=x_day)
  annual <- data.frame(YEAR=years_yr,
                       ANNUAL=x_yr)

  # compute annual sum of daily values
  daily_sum <- dplyr::group_by(daily, YEAR)
  daily_sum <- dplyr::summarise(daily_sum, DAILY_SUM=sum(VALUE))

  # compute annual ratios
  annual <- dplyr::left_join(annual, daily_sum, by='YEAR')
  annual <- dplyr::mutate(annual,
                          RATIO=ifelse(DAILY_SUM>0, ANNUAL/DAILY_SUM, 0),
                          RATIO=ifelse(RATIO>max_ratio, max_ratio, RATIO),
                          RATIO=ifelse(RATIO<min_ratio, min_ratio, RATIO))

  daily <- dplyr::left_join(daily, annual[, c('YEAR', 'RATIO')], by="YEAR")
  daily <- dplyr::mutate(daily, ADJUSTED=VALUE*RATIO)

  list(annual_ratios=annual,
       original=x_day,
       adjusted=daily[['ADJUSTED']])
}