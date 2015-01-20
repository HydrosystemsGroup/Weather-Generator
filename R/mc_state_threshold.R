#' Determine Markov State Thresholds from Precipitation
#'
#' @param x vector of daily precipitation
#' @param months vector of months corresponding to daily precipitation vector \code{x}
#' @param dry_wet_threshold threshold precipitation amount for dry/wet states
#' @param wet_extreme_quantile_threshold threshold quantile for wet/extreme states
#' @export
#' @return monthly list of precipitation thresholds as vectors of length 2 (first element is threshold amount between dry/wet, second element is threshold amount between wet/extreme)
mc_state_threshold <- function(x, months, dry_wet_threshold=0.3, wet_extreme_quantile_threshold=0.8) {
  stopifnot(any(!is.na(x)))
  stopifnot(any(!is.na(months)))
  stopifnot(length(x) == length(months))

  df <- data.frame(MONTH=months, PRCP=x)

  thresh <- plyr::dlply(df, c("MONTH"), function(x) {
    c(dry_wet_threshold, unname(quantile(x[['PRCP']], probs=wet_extreme_quantile_threshold)))
  })

  thresh
}