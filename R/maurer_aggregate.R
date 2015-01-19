#' Aggregate monthly timeseries from maurer to annual timesteps
#'
#' @param x data frame containing monthly maurer timeseries
#' @export
#' @examples
#' maurer_aggregate(maurer_get_mon(42, -71))
#'
maurer_aggregate <- function(x) {
  require(dplyr)
  x <- dplyr::mutate(x, N_DAY=lubridate::days_in_month(DATE))
  grp <- dplyr::group_by(x, LAT, LON, YEAR)
  x.yr <- dplyr::summarise(grp,
                           PRCP=sum(PRCP),
                           TMAX=sum(TMAX*N_DAY)/sum(N_DAY),
                           TMIN=sum(TMIN*N_DAY)/sum(N_DAY),
                           TAVG=sum(TAVG*N_DAY)/sum(N_DAY),
                           WIND=sum(WIND*N_DAY)/sum(N_DAY))
  x.yr <- dplyr::ungroup(x.yr)
  x.yr <- dplyr::select(x.yr, YEAR, PRCP, TMAX, TMIN, TAVG, WIND, LAT, LON)
  return(x.yr)
}