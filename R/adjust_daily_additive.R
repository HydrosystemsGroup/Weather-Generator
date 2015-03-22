#' Adjust Daily Values by Linear Addition
#'
#' @param x vector of values (e.g. daily temperature)
#' @param months vector of months corresponding to x
#' @param mean_change additive factor (scalar or vector of length 12 for monthly)
#' @export
adjust_daily_additive <- function(x, months, mean_change=0) {
  stopifnot(length(x)==length(months))

  if (length(mean_change) == 1) {
    mean_change <- rep(mean_change, 12)
  }
  stopifnot(length(mean_change)==12)

  xnew <- x
  for (m in seq(1, 12)) {
    idx <- which(months==m)
    xnew[idx] <- xnew[idx] + mean_change[m]
  }

  list(mean_change=mean_change,
       original=x,
       adjusted=xnew)
}
