#' k-Nearest Neighbors of Daily Data
#'
#' @param x annual values of climate timeseries (e.g. precipitation) as numeric array
#' @export
#' @examples
#' knn_daily(x=runif(30)*100+200, sig=0.90, noise.type="white", plot.flag=TRUE)
knn_daily <- function(cur_sim_PRCP, cur_sim_TEMP, PRCP_TODAY, TEMP_TODAY, PRCP_TOMORROW, TEMP_TOMORROW, TMAX_TOMORROW, TMIN_TOMORROW, DATE_TOMORROW, k, sd_monthly_PRCP, sd_monthly_TEMP, mean_monthly_PRCP, mean_monthly_TEMP, k1, count) {
  w_PRCP <- 100/sd_monthly_PRCP
  w_TEMP <- 10/sd_monthly_TEMP
  var_order <- 1:length(PRCP_TODAY)
  distance <- sqrt(w_PRCP*((cur_sim_PRCP-mean_monthly_PRCP) - (PRCP_TODAY-mean_monthly_PRCP))^2 + w_TEMP*((cur_sim_TEMP-mean_monthly_TEMP) - (TEMP_TODAY-mean_monthly_TEMP))^2)
  ordered_distances <- matrix(cbind(var_order,distance)[order(distance), ], ncol=2)
  K_Distances <- matrix(ordered_distances[1:k,],ncol=2)
  PROBS <- (1/row(K_Distances)[, 1]) / sum((1/row(K_Distances)[, 1]))
  selection <- sample(row(K_Distances)[,1], size=1, prob=PROBS, replace=TRUE)
  FINAL_PRCP <- PRCP_TOMORROW[K_Distances[selection,1]]
  FINAL_TEMP <- TEMP_TOMORROW[K_Distances[selection,1]]
  FINAL_TMAX <- TMAX_TOMORROW[K_Distances[selection,1]]
  FINAL_TMIN <- TMIN_TOMORROW[K_Distances[selection,1]]
  FINAL_DATE <- DATE_TOMORROW[K_Distances[selection,1]]
  return(c(FINAL_PRCP,FINAL_TEMP,FINAL_TMAX,FINAL_TMIN,FINAL_DATE))
}