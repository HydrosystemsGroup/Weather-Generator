% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arimas_simulate.R
\name{arimas_simulate}
\alias{arimas_simulate}
\title{Run simulations of multiple arima models}
\usage{
arimas_simulate(models, n, return_components = TRUE)
}
\arguments{
\item{models}{list of arima model objects (e.g. noise and wavelet components)}

\item{n}{number of simulation timesteps}

\item{return_components}{if TRUE returns list containing both the sum and indidividual components, otherwise returns just sum vector}
}
\description{
Run simulations of multiple arima models
}
\examples{
library(forecast)
models <- list(COMPONENT=auto.arima(sin(1:50)+rnorm(50)+10, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE),
               NOISE=auto.arima(rnorm(50), max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE))
arimas_simulate(models, n=40)
}

