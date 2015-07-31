#' Run Daily Weather Simulation
#'
#' @param obs_day daily historical observation dataset as \code{zoo} object with variables PRCP, TEMP, TMIN, TMAX, WIND
#' @param n_year number of simulation years
#' @param start_month initial month of the water year
#' @param start_water_year initial water year of simulation
#' @param include_leap_days include leap days in simulation time series
#' @param n_knn_annual number of years used in knn sampling algorithm
#' @param dry_wet_threshold threshold precipitation amount for dry/wet states
#' @param wet_extreme_quantile_threshold threshold quantile for wet/extreme states
#' @param adjust_annual_precip flag to adjust simulated daily precip to match simulated annual precip
#' @param annual_precip_adjust_limits range of maximum annual precip adjustment factors
#' @param dry_spell_changes adjustment factor(s) for dry spell durations (single value, or vector of length 12 for monthly)
#' @param wet_spell_changes adjustment factor(s) for wet spell durations (single value, or vector of length 12 for monthly)
#' @param prcp_mean_changes adjustment factor(s) for mean precip (single value, or vector of length 12 for monthly)
#' @param prcp_cv_changes adjustment factor(s) for precip CV (single value, or vector of length 12 for monthly)
#' @param temp_mean_changes adjustment factor(s) for mean temp (single value, or vector of length 12 for monthly)
#' @return a named list containing:
#' \item{\code{obs}}{the historical observation dataset used to train the simulation}
#' \item{\code{state_thresholds}}{monthly precipitation thresholds for defining Markov states based on the historical dataset}
#' \item{\code{transition_matrices}}{monthly transition matrices based on the historical dataset}
#' \item{\code{state_equilibria}}{monthly state equilibria probabilities}
#' \item{\code{sim}}{a data frame of the simulated daily weather}
#' @export
wgen_daily <- function(obs_day, n_year, start_month=10, start_water_year=2000, include_leap_days=FALSE,
                       n_knn_annual=100, dry_wet_threshold=0.3, wet_extreme_quantile_threshold=0.8,
                       adjust_annual_precip=TRUE, annual_precip_adjust_limits=c(0.9, 1.1),
                       dry_spell_changes=1, wet_spell_changes=1,
                       prcp_mean_changes=1, prcp_cv_changes=1, temp_mean_changes=0) {
  stopifnot(length(intersect(names(obs_day), c('PRCP', 'TEMP', 'TMIN', 'TMAX', 'WIND')))==5)

  require('zoo')
  obs_df <- as.data.frame(obs_day)
  obs_df[, 'DATE'] <- time(obs_day)
  obs_df[, 'WYEAR'] <- wyear(obs_df[['DATE']], start_month = start_month)
  obs_df[, 'MONTH'] <- month(obs_df[['DATE']])
  rownames(obs_df) <- NULL

  obs_wyr <- dplyr::group_by(obs_df, WYEAR)
  obs_wyr <- dplyr::summarize(obs_wyr,
                              PRCP=sum(PRCP),
                              TEMP=mean(TEMP),
                              TMIN=mean(TMIN),
                              TMAX=mean(TMAX),
                              WIND=mean(WIND))

  sim_annual <- sim_annual_arima(x = obs_wyr[['PRCP']], start_year = start_water_year, n_year = n_year)

  sim_day <- sim_daily_from_annual(prcp_yr = sim_annual[['out']],
                                   obs_day = zoo(x = obs_df[, c('PRCP', 'TEMP', 'TMIN', 'TMAX', 'WIND')],
                                                 order.by = obs_df[['DATE']]),
                                   obs_prcp_yr = zoo(x = obs_wyr[['PRCP']], order.by = obs_wyr[['WYEAR']]),
                                   start_month = start_month,
                                   include_leap_days = include_leap_days,
                                   n_knn_annual = n_knn_annual,
                                   dry_wet_threshold = dry_wet_threshold,
                                   wet_extreme_quantile_threshold = wet_extreme_quantile_threshold,
                                   adjust_annual_precip=adjust_annual_precip,
                                   annual_precip_adjust_limits=annual_precip_adjust_limits,
                                   dry_spell_changes=dry_spell_changes,
                                   wet_spell_changes=wet_spell_changes,
                                   prcp_mean_changes=prcp_mean_changes,
                                   prcp_cv_changes=prcp_cv_changes,
                                   temp_mean_changes=temp_mean_changes)
  sim_day
}