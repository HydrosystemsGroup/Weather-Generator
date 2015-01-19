#' Run simulation of an arima model
#'
#' @param model arima model object
#' @param n number of simulation timesteps
#' @export
#' @examples
#' library(forecast)
#' model <- auto.arima(sin(1:50)+rnorm(50)+10, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE)
#' arima_simulate(model, n=40)
arima_simulate <- function(model, n) {
  sim <- arima.sim(n=n,
                   list(ar=coef(model)[grepl('ar', names(coef(model)))],
                        ma=coef(model)[grepl('ma', names(coef(model)))]),
                   sd = sqrt(model$sigma2[[1]]))

  # extract intercept
  if ('intercept' %in% names(model$coef)) {
    intercept <- model$coef['intercept']
  } else {
    intercept <- 0
  }

  # add intercept (mean) to current simulation
  sim <- sim + intercept

  return(zoo::coredata(sim))
}


