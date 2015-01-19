#' Wavelet Analysis of Annual Timeseries
#'
#' @param x annual values of climate timeseries (e.g. precipitation) as numeric array
#' @param sig.level significance level
#' @param noise.type type of background noise ("white" or "red")
#' @export
#' @examples
#' wavelet_analysis(x=runif(30)*100+200, sig=0.90, noise.type="white")
#'
wavelet_analysis <- function(x, sig.level=0.90, noise.type=c("white", "red")) {
  noise.type <- match.arg(noise.type)
  lag1 <- switch(noise.type,
                 white = 0,
                 red = 0.72)

  bw <- biwavelet::wt(d=cbind(1:length(x), x), dt=1, dj=1/4, max.scale=length(x),
                      lag1=lag1, sig.level=sig.level, sig.test=0)

  # time-averaged global wave spectrum
  bw$gws <- apply(bw$power, 1, mean)
  bw$gws.sig <- biwavelet::wt.sig(d=x, dt=bw$dt, scale=bw$scale, sig.test=1,
                                  sig.level=sig.level, dof=length(x)-bw$scale,
                                  mother='morlet', lag1=lag1)

  return(bw)
}
