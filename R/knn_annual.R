#' Select k-Nearest Neighbor for Annual Simulation
#'
#' @param prcp target annual precipitation (scalar)
#' @param obs_prcp historical annual precipitation as zoo object
#' @param n number of years to sample
#' @export
#' @return vector of sampled years of length n
knn_annual <- function(prcp, obs_prcp, n=100) {
  stopifnot(length(prcp)==1)
  stopifnot(all(!is.na(obs_prcp)))

  k <- round(max(sqrt(length(obs_prcp)), 0.5*length(obs_prcp)), 0)
  stopifnot(k > 0)

  df <- data.frame(YEAR=time(obs_prcp),
                   PRCP=coredata(obs_prcp))

  # compute distances
  df[, 'DISTANCE'] <- sqrt((prcp - df[['PRCP']])^2)
  df <- df[order(df[['DISTANCE']]), ]

  # select k nearest
  df <- df[1:min(nrow(df), k), ]

  # compute sampling probabilities
  df[, 'ROW'] <- 1:nrow(df)
  df[, 'PROB'] <- (1/df[, 'ROW'])/(sum(1/df[, 'ROW']))
  stopifnot(abs(sum(df[['PROB']])-1) < 1e-7)

  selection <- sample(1:nrow(df), size=n, prob=df[['PROB']], replace=TRUE)
  df[['YEAR']][selection]
}
