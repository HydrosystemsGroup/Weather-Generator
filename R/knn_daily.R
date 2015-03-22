#' Select k-Nearest Neighbor for Given Simulation State
#'
#' @param x historical daily climate dataset
#' @param wday simulation water-day
#' @param state simulated state
#' @param state_prev simulated state of previous water-day
#' @param prcp_prev simulated precip of previous water-day
#' @param temp_prev simulated temp of previous water-day
#' @param prcp_sd standard deviation of historical precipitation for current month
#' @param temp_sd standard deviation of historical temperature for current month
#' @export
knn_daily <- function(x, wday, state, state_prev, prcp_prev, temp_prev, prcp_sd, temp_sd) {
  stopifnot(nrow(x)>0)

  wd_rng <- waterday_range(wday, n=7)
  current <- subset(x, WDAY %in% wd_rng & STATE==state & STATE_PREV==state_prev)

  if (nrow(current) == 0) {
    # try 30-day window
    wd_rng <- waterday_range(wday, n=30)
    current <- subset(x, WDAY %in% wd_rng & STATE==state & STATE_PREV==state_prev)
  }

  if (nrow(current) == 0) {
    # try 7-day window with any previous state
    wd_rng <- waterday_range(wday, n=7)
    current <- subset(x, WDAY %in% wd_rng & STATE==state)
  }

  if (nrow(current) == 0) {
    # try 30-day window with any previous state
    wd_rng <- waterday_range(wday, n=30)
    current <- subset(x, WDAY %in% wd_rng & STATE==state)
  }

  if (nrow(current) == 0) {
    # try full window with any previous state
    current <- subset(x, STATE==state)
  }

  if (nrow(current) == 0) {
    # usey any day in x
    current <- x
  }

  stopifnot(nrow(current)>0)
  stopifnot(all(current$STATE == state))
#   stopifnot(all(current$STATE_PREV == state_prev))

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
