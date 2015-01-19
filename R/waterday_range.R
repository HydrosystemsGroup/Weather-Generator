#' Get a range of water days around some center value
#'
#' @param x central water day
#' @param n length of range
#' @export
#' @examples
#' waterday_range(x=1, n=7)
#' waterday_range(x=10, n=7)
#' waterday_range(x=363, n=7)
waterday_range <- function(x, n=7) {
  if (n %% 2 == 1) {
    l <- (n-1)/2
    u <- (n-1)/2
  } else {
    l <- n/2
    u <- n/2 - 1
  }
  rng <- seq(x-l, x+u)
  rng <- ifelse(rng<=0, 365+rng, rng)
  rng <- ifelse(rng>365, rng-365, rng)
  rng
}
