% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wyear_date_range.R
\name{wyear_date_range}
\alias{wyear_date_range}
\title{Create Sequence of Dates by Water Year}
\usage{
wyear_date_range(n_year, start_month = 10, start_water_year = 2000,
  include_leap_days = FALSE)
}
\arguments{
\item{n_year}{number of simulation years}

\item{start_month}{first month of the water year}

\item{start_water_year}{initial water year (e.g. if start_water_year=2000 and start_month=10 then the first date is 1999-10-01)}

\item{include_leap_days}{boolean flag to include leap days (default=FALSE)}
}
\description{
Create Sequence of Dates by Water Year
}
\examples{
wyear_date_range(n_year=10, start_month=10, start_water_year=2000, include_leap_days=FALSE)
}

