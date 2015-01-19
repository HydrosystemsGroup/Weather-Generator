#' Simulate Daily Weather from Simulated Annual Precipitation
#'
#' @param prcp_yr vector of target annual precipitation
#' @param obs_da data frame of daily observed weather
#' @param obs_prcp_yr vector of observed annual precipitation
#' @param start_month first month of the water year
#' @param include_leap_days boolean flag to include leap days (default=FALSE)
#' @param n_knn_annual number of years for annual KNN sampling
#' @export
#' @return vector of sampled years of length n
#'
sim_daily_from_annual <- function(prcp_yr, obs_da, obs_prcp_yr, start_month=10, include_leap_days=FALSE, n_knn_annual=100) {
  require(zoo)
  n_year <- length(prcp_yr)
#   sim_dates <- wyear_date_range(n_year = n_year, start_month = start_month,
#                                 start_water_year = start(prcp_yr), include_leap_days = include_leap_days)
#
#   sim <- data.frame(SIM_YEAR=rep(seq(from=1, to=n_year, by=1), each=365),
#                     DATE=sim_dates,
#                     MONTH=lubridate::month(sim_dates),
#                     WDAY=rep(1:365, times=n_year))

  obs_df <- as.data.frame(obs_da)
  obs_df[, 'DATE'] <- time(obs_da)
  obs_df[, 'WYEAR'] <- wyear(obs_df[['DATE']], start_month = start_month)
  rownames(obs_df) <- NULL

  # loop through each simulation year
  sim_day_list <- lapply(seq(1, n_year), function(i) {
    prcp <- coredata(prcp_yr[i])
    yr <- time(prcp_yr[i])

    # create population of years with knn
    pop_years <- knn_annual(prcp=prcp, obs_prcp=obs_prcp_yr, n=n_knn_annual)

    # loop through population years and extract
    pop_days <- lapply(pop_years, function(yr) {
      obs_df[which(obs_df$WYEAR==yr), ]
    })

    pop_days <- do.call(rbind, pop_days)
    stopifnot(all(unique(pop_days$WYEAR) %in% unique(pop_years)))

    # perform daily markov chain simulation and knn sampling
    # to create daily sequence for current year
    sampled_days <- wgen_historical_day(pop_days, n_year = 1, start_month = start_month,
                                        start_water_year = start(prcp_yr)+i-1, include_leap_days = include_leap_days)

    # confirm that the sampled days all have the water years in the sampled years
    stopifnot(all(unique(wyear(sampled_days$SAMPLE_DATE, start_month = start_month)) %in% unique(pop_years)))

    sampled_days['SIM_YEAR'] <- i
    sampled_days
  })
  sim_day <- do.call(rbind, sim_day_list)
#   sim_day <- mutate(sim_day, DATE=DATE+years(SIM_YEAR-1))
  sim_day
}