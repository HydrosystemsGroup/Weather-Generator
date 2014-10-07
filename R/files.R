#' Save weather generator output to zip of csv files
#'
#' @param x list of results from gen_month_arma_trend()
#' @export
#' @return absolute path to zip file
#' @examples
#' clim.mon <- get_maurer_mon(lat=42, lon=-72)
#' x <- gen_month_arma_trend(x.reg.yr=clim.mon, x.loc.mon=clim.mon, n.iter=10, n.year=50, temp.factors=seq(0, 5), prcp.factors=seq(-0.5, 0.5, by=0.5))
#' save_wgen_to_zip(x)
save_wgen_to_zip <- function(x) {
  tstamp <- format(Sys.time(), '%Y%m%d_%H%M')
  temp_base <- tempdir()
  temp_path <- file.path(temp_base, tstamp)
  dir.create(temp_path)

  csv_files <- sapply(x, function(z) {
    csv_name <- paste0('climate_PRCP_', z$trend.factors$PRCP, '_TEMP_', z$trend.factors$TEMP, '.csv')
    csv_path <- file.path(temp_path, csv_name)
    write.csv(z$df, file=csv_path, row.names=FALSE)
    return(csv_path)
  })

  zip_name <- paste0('climate_', tstamp, '.zip')
  zip_path <- file.path(temp_base, zip_name)
  zip(zip_path, csv_files, flags = "-jr9X")
  return(zip_path)
}
