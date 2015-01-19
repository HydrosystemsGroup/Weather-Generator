#' Run Daily Simulation using KNN/Markov Chain
#'
#' @param n_year number of simulation years
#' @param historical historical climate dataframe
#' @param states character list of markov states
#' @param transitions list of monthly transition matrices generated with mc_fit()
#' @param start_month first month of the water year
#' @export
sim_mc_knn_day <- function(n_year, historical, states, transitions, start_month=10) {
  historical_stats <- dplyr::select(historical, MONTH, PRCP, TMAX, TMIN, TEMP)
  historical_stats <- tidyr::gather(historical_stats, VAR, VALUE, PRCP:TEMP)
  historical_stats <- plyr::dlply(historical_stats, c("VAR"), function(x) {
    plyr::dlply(x, c("MONTH"), function(x) {
      data.frame(MEAN=mean(x$VALUE), SD=sd(x$VALUE))
    })
  })

  sim <- .createSim(n_year=n_year, start_month=start_month)

  # pick random initial condition from historical dataset
  initial <- sample(which(historical$WDAY==sim[[1, 'WDAY']]), size=1)
  sim[1, 'SAMPLE_DATE'] <- historical[[initial, 'DATE']]
  sim[1, 'STATE'] <- as.character(historical[[initial, 'STATE']])
  sim[1, 'PRCP'] <- historical[[initial, 'PRCP']]
  sim[1, 'TEMP'] <- historical[[initial, 'TEMP']]
  sim[1, 'TMIN'] <- historical[[initial, 'TMIN']]
  sim[1, 'TMAX'] <- historical[[initial, 'TMAX']]

  # run markov chain simulation
  sim$STATE <- mc_sim(months=lubridate::month(sim$DATE), initial=sim[[1, 'STATE']], states=states, transition=transitions)

  for (i in 2:nrow(sim)) {
    m <- sim[i, 'MONTH']
    selected <- select_knn_day(wday=sim[[i, 'WDAY']], state=sim[[i, 'STATE']], state_prev=sim[[i-1, 'STATE']],
                               prcp_prev=sim[[i-1, 'PRCP']], temp_prev=sim[[i-1, 'TEMP']],
                               prcp_sd=historical_stats[['PRCP']][[m]][['SD']], temp_sd=historical_stats[['TEMP']][[m]][['SD']],
                               historical=historical)

    sim[i, 'SAMPLE_DATE'] <- selected[[1, 'DATE']]
    sim[i, 'PRCP'] <- selected[[1, 'PRCP']]
    sim[i, 'TEMP'] <- selected[[1, 'TEMP']]
    sim[i, 'TMIN'] <- selected[[1, 'TMIN']]
    sim[i, 'TMAX'] <- selected[[1, 'TMAX']]
  }
  sim$SAMPLE_DATE <- as.Date(sim$SAMPLE_DATE, origin='1970-01-01')

  sim
}
