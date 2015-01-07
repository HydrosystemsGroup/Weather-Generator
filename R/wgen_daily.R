.createSim <- function(n_year) {
  # TODO: leap year support? starting sim year?
  sim <- expand.grid(WDAY=seq(1, 365), WYEAR=seq(1, n_year, by=1))
  sim <- mutate(sim,
                DATE=as.Date(paste(2001, start_mon, 1, sep='-')) + lubridate::days(WDAY) - 1,
                SAMPLE_DATE=NA,
                MONTH=lubridate::month(DATE),
                STATE=NA_character_,
                PRCP=NA_real_,
                TEMP=NA_real_,
                TMIN=NA_real_,
                TMAX=NA_real_)
  sim
}

#' Run Daily Simulation using KNN/Markov Chain
#'
#' @param n_year number of simulation years
#' @param historical historical climate dataframe
#' @param states character list of markov states
#' @param transitions list of monthly transition matrices generated with mc_fit()
#' @export
sim_mc_knn_day <- function(n_year, historical, states, transitions) {
  historical_stats <- dplyr::select(historical, MONTH, PRCP, TMAX, TMIN, TEMP)
  historical_stats <- tidyr::gather(historical_stats, VAR, VALUE, PRCP:TEMP)
  historical_stats <- plyr::dlply(historical_stats, c("VAR"), function(x) {
    plyr::dlply(x, c("MONTH"), function(x) {
      data.frame(MEAN=mean(x$VALUE), SD=sd(x$VALUE))
    })
  })

  sim <- .createSim(n_year=n_year)

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

#' Select k-Nearest Neighbor for Given Simulation State
#'
#' @param wday simulation water-day
#' @param state simulated state
#' @param state_prev simulated state of previous water-day
#' @param prcp_prev simulated precip of previous water-day
#' @param temp_prev simulated temp of previous water-day
#' @param prcp_sd standard deviation of historical precipitation for current month
#' @param temp_sd standard deviation of historical temperature for current month
#' @param historical historical climate dataframe
#' @export
select_knn_day <- function(wday, state, state_prev, prcp_prev, temp_prev, prcp_sd, temp_sd, historical) {
  wd_rng <- waterday_range(wday, n=7)

  current <- subset(historical, WDAY %in% wd_rng & STATE==state & STATE_PREV==state_prev)
  stopifnot(nrow(current)>1)
  stopifnot(all(current$STATE == state))
  stopifnot(all(current$STATE_PREV == state_prev))

  k <- round(sqrt(nrow(current)))

  w_PRCP <- 100/prcp_sd
  w_TEMP <- 10/temp_sd

  current$DISTANCE <- sqrt(w_PRCP*(prcp_prev - current[['PRCP_PREV']])^2 + w_TEMP*(temp_prev - current[['TEMP_PREV']])^2)

  current <- dplyr::arrange(current, DISTANCE)
  current <- dplyr::filter(current, dplyr::row_number(DISTANCE)<=k)
  current <- dplyr::mutate(current, RANK=dplyr::row_number(DISTANCE), PROB=(1/RANK)/sum(1/RANK))

  selection <- sample(1:nrow(current), size=1, prob=current$PROB)
  current[selection, ]
}


#' Run Daily Weather Generator using Historical Data
#'
#' @param historical historical climate dataframe
#' @param n_year number of simulation years
#' @param start_mon initial month of the water year
#' @export
wgen_historical_day <- function(historical, n_year, start_mon=10) {
  # add MONTH and WDAY columns
  historical <- dplyr::mutate(historical, MONTH=month(DATE), WDAY=waterday(DATE, start.month=start_mon))

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