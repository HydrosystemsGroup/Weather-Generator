% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arima_simulate.R
\name{arima_simulate}
\alias{arima_simulate}
\title{Run simulation of an arima model}
\usage{
arima_simulate(model, n)
}
\arguments{
\item{model}{arima model object}

\item{n}{number of simulation timesteps}
}
\description{
Run simulation of an arima model
}
\examples{
library(forecast)
model <- auto.arima(sin(1:50)+rnorm(50)+10, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE)
arima_simulate(model, n=40)
}

