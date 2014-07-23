#' k-Nearest Neighbors of monthly climate
#'
#' @param cur_sim_PRCP value of simulated daily precipitation
#' @param cur_sim_TEMP value of simulated daily mean air temperature
#' @param PRCP_TODAY amount of precipitation for today
#' @param TEMP_TODAY mean air temperature for today
#' @param PRCP_TOMORROW amount of precipitation for tomorrow
#' @param TEMP_TOMORROW mean air temperature for tomorrow
#' @param TMAX_TOMORROW maximum air temperature for tomorrow
#' @param TMIN_TOMORROW minimum air temperature for tomorrow
#' @param DATE_TOMORROW date of tomorrow
#' @param k number of nearest neighbors
#' @param sd_monthly_PRCP standard deviation of monthly precip
#' @param sd_monthly_TEMP standard deviation of monthly air temperature
#' @param mean_monthly_PRCP mean of monthly precip
#' @param mean_monthly_TEMP mean of monthly air temperature
#' @param k1 trial number
#' @param count current timestep
#' @export
#' @examples
#' KNN(cur_sim_PRCP,cur_sim_TEMP,PRCP_TODAY,TEMP_TODAY,PRCP_TOMORROW,TEMP_TOMORROW,TMAX_TOMORROW,TMIN_TOMORROW,DATE_TOMORROW,k,sd_monthly_PRCP,sd_monthly_TEMP,mean_monthly_PRCP,mean_monthly_TEMP,k1,count)

KNN <- function(cur_sim_PRCP,cur_sim_TEMP,PRCP_TODAY,TEMP_TODAY,PRCP_TOMORROW,TEMP_TOMORROW,TMAX_TOMORROW,TMIN_TOMORROW,DATE_TOMORROW,k,sd_monthly_PRCP,sd_monthly_TEMP,mean_monthly_PRCP,mean_monthly_TEMP,k1,count) {
  w_PRCP <- 100/sd_monthly_PRCP
  w_TEMP <- 10/sd_monthly_TEMP
  var_order <- 1:length(PRCP_TODAY)
  distance <- sqrt(w_PRCP*((cur_sim_PRCP-mean_monthly_PRCP) - (PRCP_TODAY-mean_monthly_PRCP))^2 + w_TEMP*((cur_sim_TEMP-mean_monthly_TEMP) - (TEMP_TODAY-mean_monthly_TEMP))^2)
  ordered_distances <- matrix(cbind(var_order,distance)[order(distance),],ncol=2)
  K_Distances <- matrix(ordered_distances[1:k,],ncol=2)
  PROBS <- (1/row(K_Distances)[,1]) / sum((1/row(K_Distances)[,1]))
  set.seed(k1*count)
  selection <- sample(row(K_Distances)[,1],size=1,prob=PROBS,replace=TRUE)
  FINAL_PRCP <- PRCP_TOMORROW[K_Distances[selection,1]]
  FINAL_TEMP <- TEMP_TOMORROW[K_Distances[selection,1]]
  FINAL_TMAX <- TMAX_TOMORROW[K_Distances[selection,1]]
  FINAL_TMIN <- TMIN_TOMORROW[K_Distances[selection,1]]
  FINAL_DATE <- DATE_TOMORROW[K_Distances[selection,1]]
  return(c(FINAL_PRCP,FINAL_TEMP,FINAL_TMAX,FINAL_TMIN,FINAL_DATE))
}
