#' Determine Markov State Thresholds from Precipitation
#'
#' @param prcp vector of daily precipitation
#' @param months vector of months corresponding to daily precipitation
#' @param dry_wet threshold precipitation amount for dry/wet states
#' @param wet_extreme_quantile threshold quantile for wet/extreme states
#' @export
mc_state_threshold <- function(prcp, months, dry_wet=0.3, wet_extreme_quantile=0.8) {
  stopifnot(any(!is.na(prcp)))
  stopifnot(any(!is.na(months)))
  stopifnot(length(prcp) == length(months))

  df <- data.frame(MONTH=months, PRCP=prcp)

  thresh <- plyr::dlply(df, c("MONTH"), function(x) {
    c(dry_wet, unname(quantile(x[['PRCP']], probs=wet_extreme_quantile)))
  })

  thresh
}