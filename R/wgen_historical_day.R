#' Run Daily Weather Generator using Historical Data
#'
#' @param historical historical climate dataframe
#' @param n_year number of simulation years
#' @param start_month initial month of the water year
#' @export
wgen_historical_day <- function(historical, n_year, start_month=10) {
  # add MONTH and WDAY columns
  historical <- dplyr::mutate(historical, MONTH=month(DATE), WDAY=waterday(DATE, start.month=start_month))

  # compute precipitation thresholds by month
  thresh <- mc_state_threshold(historical[['PRCP']], historical[['MONTH']],
                               dry_wet=0.3, wet_extreme_quantile=0.8)

  states <- c('d', 'w', 'e')

  # assign precipitation state and lagged variable columns
  historical <- dplyr::mutate(historical,
                              STATE=mc_assign_states(PRCP, MONTH, states, thresh),
                              STATE_PREV=lag(STATE),
                              PRCP_PREV=lag(PRCP),
                              TEMP_PREV=lag(TEMP),
                              TMAX_PREV=lag(TMAX),
                              TMIN_PREV=lag(TMIN))

  # fit Markov Chain transition probabilities
  transitions <- mc_fit(states=historical[['STATE']], months=historical[['MONTH']])

  # run daily simulation
  sim <- sim_mc_knn_day(n_year=n_year, historical=historical, states=states, transitions=transitions)

  sim
}
