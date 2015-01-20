#' Run Daily Weather Simulation
#'
#' @param historical historical climate dataframe
#' @param n_year number of simulation years
#' @param dry_wet_threshold threshold precipitation amount for dry/wet states
#' @param wet_extreme_quantile_threshold threshold quantile for wet/extreme states
#' @param start_month initial month of the water year
#' @param start_water_year initial water year of simulation
#' @param include_leap_days include leap days in simulation time series
#' @return a named list containing:
#' \item{\code{x}}{the historical dataset used to train the simulation}
#' \item{\code{state_thresholds}}{monthly precipitation thresholds for defining Markov states based on the historical dataset}
#' \item{\code{transition_matrices}}{monthly transition matrices based on the historical dataset}
#' \item{\code{state_equilibria}}{monthly state equilibria probabilities}
#' \item{\code{out}}{a data frame of the simulated daily weather}
#' @export
sim_daily <- function(historical, n_year, dry_wet_threshold=0.3, wet_extreme_quantile_threshold=0.8,
                      start_month=10, start_water_year=2000, include_leap_days=FALSE,
                      dry_spell_changes=1, wet_spell_changes=1) {

  if (length(dry_spell_changes) == 1) {
    dry_spell_changes <- rep(dry_spell_changes, 12)
  }
  if (length(wet_spell_changes) == 1) {
    wet_spell_changes <- rep(wet_spell_changes, 12)
  }
  stopifnot(length(dry_spell_changes)==12)
  stopifnot(length(wet_spell_changes)==12)

  # add MONTH and WDAY columns
  historical <- dplyr::mutate(historical, MONTH=month(DATE), WDAY=waterday(DATE, start_month=start_month))

  # compute precipitation thresholds by month
  states <- c('d', 'w', 'e')
  thresh <- mc_state_threshold(historical[['PRCP']], historical[['MONTH']],
                               dry_wet_threshold=0.3, wet_extreme_quantile_threshold=0.8)

  # assign precipitation state and lagged variable columns
  historical <- dplyr::mutate(historical,
                              STATE=as.character(mc_assign_states(PRCP, MONTH, states, thresh)),
                              STATE_PREV=lag(STATE),
                              PRCP_PREV=lag(PRCP),
                              TEMP_PREV=lag(TEMP),
                              TMAX_PREV=lag(TMAX),
                              TMIN_PREV=lag(TMIN))

  # fit Markov Chain transition probabilities
  transitions_historical <- mc_fit(states=historical[['STATE']], months=historical[['MONTH']])
  equilibria_historical <- lapply(transitions_historical, mc_state_equilibrium)
  probability_wet_historical <- sapply(equilibria_historical, function(equil) {
    sum(equil[2:3])
  })

  # adjust Markov Chain transition probabilities to change dry/wet spell
  transitions <- lapply(1:length(transitions_historical), function(i) {
    x <- mc_adjust_transition(transitions_historical[[i]],
                              dry_spell=dry_spell_changes[[i]],
                              wet_spell=wet_spell_changes[[i]])
    x
  })
  equilibria <- lapply(transitions, mc_state_equilibrium)
  probability_wet <- sapply(equilibria, function(equil) {
    sum(equil[2:3])
  })

  ratio_probability_wet <- probability_wet_historical/probability_wet

  # run daily simulation
  sim <- sim_mc_knn_day(n_year=n_year, historical=historical, states=states, transitions=transitions,
                        start_month=start_month, start_water_year=start_water_year,
                        include_leap_days=include_leap_days)
  sim_prcp_sampled <- select(sim, DATE, PRCP)

  list(x=historical,
       state_thresholds=thresh,
       transition_matrices=transitions,
       state_equilibria=equilibria,
       ratio_probability_wet=ratio_probability_wet,
       out=sim)
}
