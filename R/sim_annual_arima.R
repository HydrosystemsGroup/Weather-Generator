#' Simulate Annual Timeseries with ARIMA Model
#'
#' Fit an ARIMA model to a given vector of annual precipitation \code{x}, and run a simulation of length \code{n_year} beginning with year \code{start_year}.
#'
#' @param x vector of observed annual values
#' @param start_year first simulation year
#' @param n_year number of simulation years
#' @export
#' @return a named list containing:
#' \item{\code{model}}{the fitted arima model returned by \code{forecast::auto.arima()}}
#' \item{\code{x}}{the historical annual precipitation as a numeric vector}
#' \item{\code{out}}{the simulated annual precipitation of length \code{n_year} as a \code{zoo} object}
#' @examples
#' sim_annual_arima(rnorm(100))
sim_annual_arima <- function(x, start_year=2000, n_year=10) {
  stopifnot(all(!is.na(x)))

  # fit ARIMA model
  ar_model <- forecast::auto.arima(x, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE)

  # simulation ARIMA model
  sim_x <- arima_simulate(model=ar_model, n=n_year)

  # create output zoo object
  sim_years <- seq(start_year, by = 1, length.out = n_year)
  out <- zoo::zoo(x = sim_x, order.by = sim_years)

  # return list
  list(model=ar_model,
       x=x,
       out=out)
}
