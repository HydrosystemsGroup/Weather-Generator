#' Run Monte Carlo Simulation of Multiple ARIMA Models
#'
#' @param models list of arima model objects (e.g. noise and wavelet components)
#' @param n number of simulation timesteps
#' @param n.iter number of monte carlo iterations
#' @export
#' @examples
#' library(forecast)
#' models <- list(COMPONENT=auto.arima(sin(1:50)+rnorm(50)+10, max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE),
#'                NOISE=auto.arima(rnorm(50), max.p=2, max.q=2, max.P=0, max.Q=0, stationary=TRUE))
#' arimas_monte(models, n=40, n.iter=100)
arimas_monte <- function(models, n, n.iter) {
  # run dummy simulation and wavelet analysis to get number of periods for wavelet transform
  temp.sim.warm <- simulate_arimas(models=models, n=n)
  temp.wt <- wavelet_analysis(temp.sim.warm$sum, sig.level=0.90, noise.type='white')

  sim.x <- array(NA, c(n, n.iter))
  sim.gws <- array(NA, c(length(temp.wt$gws), n.iter))
  for (i in 1:n.iter) {
    i.sim.warm <- simulate_arimas(models=models, n=n)
    i.wt <- wavelet_analysis(i.sim.warm$sum, sig.level=0.90, noise.type='white')
    sim.x[,i] <- i.sim.warm$sum
    sim.gws[,i] <- i.wt$gws
  }

  sim.x.stat <- cbind(MEAN=rowMeans(sim.x),
                      Q025=apply(sim.x, 1, quantile, 0.025),
                      Q975=apply(sim.x, 1, quantile, 0.975))

  sim.gws.stat <- cbind(MEAN=rowMeans(sim.gws),
                        Q025=apply(sim.gws, 1, quantile, 0.025),
                        Q975=apply(sim.gws, 1, quantile, 0.975))

  return(list(x=sim.x, gws=sim.gws, x.stat=sim.x.stat, gws.stat=sim.gws.stat))
}