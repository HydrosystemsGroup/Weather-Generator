#' Daily Weather Generator Conditioned on Simulated Annual Precipitation
#'
#' Given a timeseries of simulated annual precipitation, simulate daily weather using the KNN and Markov Chain method.
#'
#' @param prcp_yr simulated annual precipitation as \code{zoo} object
#' @param obs_day observed daily weather as multi-column \code{zoo} object
#' @param obs_prcp_yr observed annual precipitation as numeric vector
#' @param start_month first month of the water year
#' @param include_leap_days boolean flag to include leap days (default=FALSE)
#' @param n_knn_annual number of years used for annual KNN sampling
#' @param dry_wet_threshold threshold precipitation amount for dry/wet states
#' @param wet_extreme_quantile_threshold threshold quantile for wet/extreme states
#' @export
#' @return a named list containing:
#' \item{\code{x}}{the historical dataset used to train the simulation}
#' \item{\code{state_thresholds}}{monthly precipitation thresholds for defining Markov states based on the historical dataset}
#' \item{\code{transition_matrices}}{monthly transition matrices based on the historical dataset}
#' \item{\code{state_equilibria}}{monthly state equilibria probabilities}
#' \item{\code{out}}{a data frame of the simulated daily weather}
#'
sim_daily_from_annual <- function(prcp_yr, obs_day, obs_prcp_yr, start_month=10, include_leap_days=FALSE,
                                  n_knn_annual=100, dry_wet_threshold=0.3, wet_extreme_quantile_threshold=0.8,
                                  adjust_annual_precip=TRUE, annual_precip_adjust_limits=c(0.9, 1.1),
                                  dry_spell_changes=1, wet_spell_changes=1,
                                  prcp_mean_changes=1, prcp_cv_changes=1, temp_changes=0) {
  if (include_leap_days) {
    stop('Leap days not currently supported')
  }

  if (length(prcp_mean_changes) == 1) {
    prcp_mean_changes <- rep(prcp_mean_changes, 12)
  }
  if (length(prcp_cv_changes) == 1) {
    prcp_cv_changes <- rep(prcp_cv_changes, 12)
  }
  if (length(temp_changes) == 1) {
    temp_changes <- rep(temp_changes, 12)
  }
  stopifnot(length(prcp_mean_changes)==12)
  stopifnot(length(prcp_cv_changes)==12)
  stopifnot(length(temp_changes)==12)

  n_year <- length(prcp_yr)

  require('zoo')
  obs_df <- as.data.frame(obs_day)
  obs_df[, 'DATE'] <- time(obs_day)
  obs_df[, 'WYEAR'] <- wyear(obs_df[['DATE']], start_month = start_month)
  obs_df[, 'MONTH'] <- month(obs_df[['DATE']])
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
    sim_days_i <- sim_daily(pop_days, n_year = 1, start_month = start_month,
                            start_water_year = start(prcp_yr)+i-1, include_leap_days = include_leap_days,
                            dry_spell_changes = dry_spell_changes,
                            wet_spell_changes = wet_spell_changes)

    # confirm that the sampled days all have the water years in the sampled years
    stopifnot(all(unique(wyear(sim_days_i[['out']][['SAMPLE_DATE']], start_month = start_month)) %in% unique(pop_years)))

    sim_days_i[['out']]['SIM_YEAR'] <- i
    sim_days_i
  })

  sim_day_list_out <- lapply(sim_day_list, function(x) {
    x[['out']]
  })

  ratio_probability_wet <- lapply(sim_day_list, function(x) {
    x[['ratio_probability_wet']]
  })
  ratio_probability_wet <- do.call(cbind, ratio_probability_wet)
  mean_ratio_probability_wet <- apply(ratio_probability_wet, 1, mean)
  stopifnot(length(mean_ratio_probability_wet)==12)

  sim_day <- do.call(rbind, sim_day_list_out)

  sum_unadjusted <- sim_day
  if (adjust_annual_precip) {
    annual_sum_adjustment <- adjust_daily_to_annual(values_day=sim_day[['PRCP']],
                                                    years_day=wyear(sim_day[['DATE']], start_month = start_month),
                                                    values_yr=zoo::coredata(prcp_yr),
                                                    years_yr=time(prcp_yr),
                                                    min_ratio=annual_precip_adjust_limits[1],
                                                    max_ratio=annual_precip_adjust_limits[2])
    sim_day[, 'PRCP'] <- annual_sum_adjustment$adjusted
  }

  gamma_adjustment <- adjust_daily_gamma(sim_day[['PRCP']], sim_day[['MONTH']],
                                         mean_change=mean_ratio_probability_wet*prcp_mean_changes,
                                         cv_change=prcp_cv_changes)
  sim_day[, 'PRCP'] <- gamma_adjustment$adjusted

  tmin_adjustment <- adjust_daily_additive(x = sim_day[['TMIN']],
                                           months = sim_day[['MONTH']],
                                           changes = temp_changes)
  tmax_adjustment <- adjust_daily_additive(x = sim_day[['TMAX']],
                                           months = sim_day[['MONTH']],
                                           changes = temp_changes)
  temp_adjustment <- adjust_daily_additive(x = sim_day[['TEMP']],
                                           months = sim_day[['MONTH']],
                                           changes = temp_changes)
  sim_day[, 'TMIN'] <- tmin_adjustment$adjusted
  sim_day[, 'TMAX'] <- tmax_adjustment$adjusted
  sim_day[, 'TEMP'] <- temp_adjustment$adjusted

  # compute precipitation thresholds by month
  thresh <- mc_state_threshold(obs_df[['PRCP']], obs_df[['MONTH']],
                               dry_wet_threshold=dry_wet_threshold,
                               wet_extreme_quantile_threshold=wet_extreme_quantile_threshold)
  obs_df[, 'STATE'] <- mc_assign_states(obs_df[['PRCP']], obs_df[['MONTH']], states=c('d', 'w', 'e'), thresh=thresh)

  # compute monthly transition matrices
  transitions <- mc_fit(states=obs_df[['STATE']], months=obs_df[['MONTH']])
  equilibria <- lapply(transitions, mc_state_equilibrium)

  list(x=obs_df,
       state_thresholds=thresh,
       transition_matrices=transitions,
       state_equilibria=equilibria,
       changes=list(ratio_probability_wet=mean_ratio_probability_wet,
                    prcp_mean=prcp_mean_changes,
                    prcp_cv=prcp_cv_changes,
                    temp=temp_changes),
       out=sim_day)
}