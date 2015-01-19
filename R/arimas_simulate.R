#' Run simulations of multiple arima models
#'
#' @param models list of arima model objects (e.g. noise and wavelet components)
#' @param n number of simulation timesteps
#' @param return_components if TRUE returns list containing both the sum and indidividual components, otherwise returns just sum vector
#' @export
#' @examples
#' library(forecast)
#' models <- list(COMPONENT=auto.arima(sin(1:50)+rnorm(50)+10, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE),
#'                NOISE=auto.arima(rnorm(50), max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE))
#' arimas_simulate(models, n=40)
arimas_simulate <- function(models, n, return_components=TRUE) {
  # run simulation on each model and combine into 2-d array
  sim.components <- lapply(models, arima_simulate, n=n)
  sim.components <- sapply(sim.components, cbind)

  # compute sum of individual component simulations
  if (ncol(sim.components) > 1) {
    sim.sum <- apply(sim.components, 1, FUN=sum)
  } else {
    sim.sum <- sim.components
  }

  if (return_components) {
    return(list(components=sim.components, sum=sim.sum))
  } else {
    return(sim.sum)
  }
}
