#' Run Daily Simulation using KNN/Markov Chain
#'
#' @param x historical daily climate dataset
#' @param n_year number of simulation years
#' @param states character list of markov states
#' @param transitions list of monthly transition matrices generated with mc_fit()
#' @param start_month first month of the water year
#' @param start_water_year initial water year of simulation
#' @param include_leap_days include leap days in simulation time series
#' @export
sim_mc_knn_day <- function(x, n_year, states, transitions,
                           start_month=10, start_water_year=2000, include_leap_days=FALSE) {
  x_stats <- dplyr::select(x, MONTH, PRCP, TMAX, TMIN, TEMP)
  x_stats <- tidyr::gather(x_stats, VAR, VALUE, PRCP:TEMP)
  x_stats <- plyr::dlply(x_stats, c("VAR"), function(x) {
    plyr::dlply(x, c("MONTH"), function(x) {
      data.frame(MEAN=mean(x$VALUE), SD=sd(x$VALUE))
    })
  })

  sim_dates <- wyear_date_range(n_year = n_year, start_month = start_month,
                                start_water_year = start_water_year, include_leap_days = include_leap_days)
  sim <- data.frame(SIM_YEAR=rep(seq(from=1, to=n_year, by=1), each=365),
                    DATE=sim_dates,
                    MONTH=lubridate::month(sim_dates),
                    WDAY=rep(1:365, times=n_year))

  # pick random initial condition from x dataset
  initial <- sample(which(x$WDAY==sim[[1, 'WDAY']]), size=1)
  sim[1, 'SAMPLE_DATE'] <- format(x[[initial, 'DATE']], "%Y-%m-%d")
  sim[1, 'STATE'] <- as.character(x[[initial, 'STATE']])
  sim[1, 'PRCP'] <- x[[initial, 'PRCP']]
  sim[1, 'TEMP'] <- x[[initial, 'TEMP']]
  sim[1, 'TMIN'] <- x[[initial, 'TMIN']]
  sim[1, 'TMAX'] <- x[[initial, 'TMAX']]
  sim[1, 'WIND'] <- x[[initial, 'WIND']]

  # run markov chain simulation
  sim$STATE <- mc_simulate(months=lubridate::month(sim$DATE), initial=sim[[1, 'STATE']], states=states, transition=transitions)

  for (i in 2:nrow(sim)) {
    m <- sim[i, 'MONTH']
    selected <- knn_daily(x=x, wday=sim[[i, 'WDAY']],
                          state=sim[[i, 'STATE']], state_prev=sim[[i-1, 'STATE']],
                          prcp_prev=sim[[i-1, 'PRCP']], temp_prev=sim[[i-1, 'TEMP']],
                          prcp_sd=x_stats[['PRCP']][[m]][['SD']], temp_sd=x_stats[['TEMP']][[m]][['SD']])

    sim[i, 'SAMPLE_DATE'] <- format(selected[[1, 'DATE']], "%Y-%m-%d")
    sim[i, 'PRCP'] <- selected[[1, 'PRCP']]
    sim[i, 'TEMP'] <- selected[[1, 'TEMP']]
    sim[i, 'TMIN'] <- selected[[1, 'TMIN']]
    sim[i, 'TMAX'] <- selected[[1, 'TMAX']]
    sim[i, 'WIND'] <- selected[[1, 'WIND']]
  }
  sim$SAMPLE_DATE <- as.Date(sim$SAMPLE_DATE)

  sim
}
