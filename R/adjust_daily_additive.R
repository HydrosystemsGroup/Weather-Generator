#' Adjust Daily Values by Linear Addition
#'
#' @param x vector of values (e.g. daily temperature)
#' @param months vector of months corresponding to x
#' @param changes additive factor (scalar or vector of length 12 for monthly)
#' @export
adjust_daily_additive <- function(x, months, changes=0) {
  stopifnot(length(x)==length(months))

  if (length(changes) == 1) {
    changes <- rep(changes, 12)
  }
  stopifnot(length(changes)==12)

  xnew <- x
  for (m in seq(1, 12)) {
    idx <- which(months==m)
    xnew[idx] <- xnew[idx] + changes[m]
  }

  list(changes=changes,
       original=x,
       adjusted=xnew)
}
