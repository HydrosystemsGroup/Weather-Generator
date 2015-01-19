#' Assign Markov States from Precipitation and State Thresholds
#'
#' @param prcp vector of daily precipitation
#' @param months vector of months corresponding to daily precipitation
#' @param states character list of markov states
#' @param thresholds list of monthly transition matrices generated from mc_fit()
#' @export
mc_assign_states <- function(prcp, months, states, thresholds) {
  stopifnot(length(prcp)==length(months))
  stopifnot(length(thresholds)==12)

  x <- data.frame(PRCP=prcp, MONTH=months)
  x <- group_by(x, MONTH)
  x <- mutate(x, STATE=cut(PRCP, breaks=c(0, thresholds[[unique(MONTH)]], Inf), include.lowest=TRUE, right=TRUE, labels=states))
  x <- ungroup(x)
  x[['STATE']] <- ordered(x[['STATE']], levels=states)
  x[['STATE']]
}
