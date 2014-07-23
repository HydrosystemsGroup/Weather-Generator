#' k-Nearest Neighbors of Annual Sequence
#'
#' @param sim_annual_prcp value of simulated annual precipitation
#' @param ANNUAL_PRCP numeric array of observed annual precipitation
#' @param WATER_YEAR_A integer array of water years corresponding to observed annual precipitation
#' @param kk number of nearest neighbors (?)
#' @param k1 trial number as integer
#' @param y integer representing index of year simulation in 1:num_years_sim
#' @export
#' @examples
#' KNN_ANNUAL(sim_annual_prcp,ANNUAL_PRCP,WATER_YEAR_A,kk,k1,y)
KNN_ANNUAL <- function(sim_annual_prcp,ANNUAL_PRCP,WATER_YEAR_A,kk,k1,y) {
  var_order <- 1:length(ANNUAL_PRCP)
  distance <- sqrt((sim_annual_prcp - ANNUAL_PRCP)^2)
  ordered_distances <- matrix(cbind(var_order,distance)[order(distance),],ncol=2)
  K_Distances <- matrix(ordered_distances[1:kk,],ncol=2)
  PROBS <- (1/row(K_Distances)[,1]) / sum((1/row(K_Distances)[,1]))
  set.seed(k1*y)
  selection <- sample(row(K_Distances)[,1],size=100,prob=PROBS,replace=TRUE)
  FINAL_YEARS <- WATER_YEAR_A[K_Distances[selection,1]]
}
