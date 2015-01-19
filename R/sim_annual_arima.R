#' Simulate Annual Timeseries with ARIMA Model
#'
#' @param x vector of observed annual values
#' @param start_water_year first simulation water year
#' @param n_year number of simulation years
#' @export
#' @return list containing ARIMA model, observed values, and simulated values as zoo object
sim_annual_arima <- function(x, start_water_year=2000, n_year=10) {
  stopifnot(all(!is.na(x)))

  ar_model <- forecast::auto.arima(x, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE)
  sim_prcp <- arima_simulate(model=ar_model, n=n_year)
  sim_wyr <- seq(start_water_year, by = 1, length.out = n_year)
  sim <- zoo::zoo(x = sim_prcp, order.by = sim_wyr)

  list(model=ar_model,
       obs=x,
       sim=sim)
}
