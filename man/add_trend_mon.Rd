% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_month.R
\name{add_trend_mon}
\alias{add_trend_mon}
\title{Add climate trends to a single generated monthly climate timeseries}
\usage{
add_trend_mon(sim.gen, temp.factor, prcp.factor)
}
\arguments{
\item{sim.gen}{data.frame of a single monthly climate timeseries generated by gen_month_arma()}

\item{temp.factor}{additive temperature trend factor}

\item{prcp.factor}{multiplicative precipitation trend factor}
}
\description{
Add climate trends to a single generated monthly climate timeseries
}
\examples{
sim.gen <- gen_month_arma(x.reg.yr=clim.reg.yr, x.loc.mon=clim.loc.mon, n.iter=50, n.year=nrow(clim.reg.yr))
add_trend_mon(sim.gen=sim.gen$sim.mon[[1]], temp.factor=2, prcp.factor=1.5)

}

