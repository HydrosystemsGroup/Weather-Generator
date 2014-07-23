#' Run daily weather generator
#'
#' @param k1 trial number as integer
#' @param num_year_sim number of years for simulation
#' @param PRCP_FINAL_ANNUAL_SIM numeric array of simulated annual precipitation
#' @param ANNUALPRCP numeric array of observed annual precipitation
#' @param WATER_YEAR_A numeric array of water years for annual observed data
#' @param WATER_YEAR_D array of Dates corresponding to Jan 1 of each water year in WATER_YEAR_A
#' @param PRCP numeric array of daily precipitation (mm/day?)
#' @param TEMP numeric array of daily mean air temperature (degC?)
#' @param TMAX numeric array of daily maximum air temperature (degC?)
#' @param TMIN numeric array of daily minimum air temperature (degC?)
#' @param DATE_D array of Dates corresponding to daily input data (same length as PRCP, TEMP, TMAX, TMIN)
#' @param MONTH_D numeric array of months corresponding to DATE_D
#' @param YEAR_D numeric array of years corresponding to DATE_D
#' @param MONTH_DAY_D numeric matrix of month and day correspond to DATE_D with column names (MONTH_D, DAY_D)
#' @param wet_spell_change numeric array of monthly wet spell adjustment factors
#' @param dry_spell_change numeric array of monthly dry spell adjustment factors
#' @param month_list integer array of month sequence (starts with 10 for water year, 1 for calendar year)
#' @param water_year_start initial water year of observed data as integer
#' @param water_year_end final water year of observed data as integer
#' @export
#' @examples
#' DAILY_WEATHER_GENERATOR(k1,num_year_sim,PRCP_FINAL_ANNUAL_SIM,ANNUAL_PRCP,WATER_YEAR_A,WATER_YEAR_D,PRCP,TEMP,TMAX,TMIN,DATE_D,MONTH_D,YEAR_D,MONTH_DAY_D,wet_spell_change,dry_spell_change,month_list,water_year_start,water_year_end)
DAILY_WEATHER_GENERATOR <- function(k1,num_year_sim,PRCP_FINAL_ANNUAL_SIM,ANNUAL_PRCP,WATER_YEAR_A,WATER_YEAR_D,PRCP,TEMP,TMAX,TMIN,DATE_D,MONTH_D,YEAR_D,MONTH_DAY_D,wet_spell_change,dry_spell_change,month_list,water_year_start,water_year_end) {
  START_YEAR_SIM <- 2000
  END_YEAR_SIM <- START_YEAR_SIM + num_year_sim
  DATE_SIM <- seq(as.Date(paste(START_YEAR_SIM,"-1-01",sep="")),as.Date(paste(END_YEAR_SIM,"-12-31",sep="")),by="day")
  DAY_SIM <- as.numeric(format(DATE_SIM,"%d"))
  MONTH_SIM <- as.numeric(format(DATE_SIM,"%m"))
  YEAR_SIM <- as.numeric(format(DATE_SIM,"%Y"))
  no_leap <- which(MONTH_SIM!=2 | DAY_SIM!= 29)
  DAY_SIM <- DAY_SIM[no_leap]
  MONTH_SIM <- MONTH_SIM[no_leap]
  YEAR_SIM <- YEAR_SIM[no_leap]
  WATER_YEAR_SIM <- YEAR_SIM
  if (water_yr_change) {
    WATER_YEAR_SIM[which(MONTH_SIM>=month_list[1])] <- WATER_YEAR_SIM[which(MONTH_SIM>=month_list[1])] + 1
  }
  WATER_YEAR_LOCATIONS_SIM <- which(WATER_YEAR_SIM>=(START_YEAR_SIM+1) & WATER_YEAR_SIM<=END_YEAR_SIM)
  WATER_YEAR_SIM <- WATER_YEAR_SIM[WATER_YEAR_LOCATIONS_SIM]
  YEAR_SIM <- YEAR_SIM[WATER_YEAR_LOCATIONS_SIM]
  MONTH_SIM <- MONTH_SIM[WATER_YEAR_LOCATIONS_SIM]
  DAY_SIM <- DAY_SIM[WATER_YEAR_LOCATIONS_SIM]
  DATE_SIM <- as.Date(paste(WATER_YEAR_SIM,MONTH_SIM,DAY_SIM,sep="-"))
  SIM_LENGTH <- length(DATE_SIM)
  DATE_M_SIM <- subset(DATE_SIM,DAY_SIM==1)
  YEAR_M_SIM <- as.numeric(format(DATE_M_SIM,"%Y"))
  MONTH_M_SIM <- as.numeric(format(DATE_M_SIM,"%m"))
  DATE_A_SIM <- subset(DATE_M_SIM,MONTH_M_SIM==month_list[1])
  WATER_YEAR_A_SIM <- as.numeric(format(DATE_A_SIM,"%Y"))

  p00_final <- array(NA,SIM_LENGTH)
  p01_final <- array(NA,SIM_LENGTH)
  p02_final <- array(NA,SIM_LENGTH)
  p10_final <- array(NA,SIM_LENGTH)
  p11_final <- array(NA,SIM_LENGTH)
  p12_final <- array(NA,SIM_LENGTH)
  p20_final <- array(NA,SIM_LENGTH)
  p21_final <- array(NA,SIM_LENGTH)
  p22_final <- array(NA,SIM_LENGTH)
  OCCURENCES <- array(0,c(SIM_LENGTH))
  SIM_PRCP <- array(0,c(SIM_LENGTH))
  SIM_TEMP <- array(25,c(SIM_LENGTH))
  SIM_TMAX <- array(30,c(SIM_LENGTH))
  SIM_TMIN <- array(20,c(SIM_LENGTH))
  SIM_DATE <- array(as.Date(paste(sample(WATER_YEAR_A,size=1),month_list[1],"01",sep="-")),c(SIM_LENGTH))

  Ratio_PI_old_PI_new <- array(NA,c(num_year_sim,12))

  count <- 1
  set.seed(k1)
  rn_all <- runif(SIM_LENGTH,0,1)
  kk <- max(round(sqrt(length(ANNUAL_PRCP)),0),round(length(ANNUAL_PRCP),0)*.5)
  for (y in 1:num_year_sim) {

    sim_annual_prcp <- PRCP_FINAL_ANNUAL_SIM[y]
    CUR_YEARS <- KNN_ANNUAL(sim_annual_prcp,ANNUAL_PRCP,WATER_YEAR_A,kk,k1,y)
    conditional_selection <- NULL
    for (yy in 1:length(CUR_YEARS)) {
      conditional_selection <- c(conditional_selection,which(WATER_YEAR_D==CUR_YEARS[yy]))
    }
    PRCP_CURRENT <- PRCP[conditional_selection]
    TEMP_CURRENT <- TEMP[conditional_selection]
    TMAX_CURRENT <- TMAX[conditional_selection]
    TMIN_CURRENT <- TMIN[conditional_selection]
    DATE_D_CURRENT <- DATE_D[conditional_selection]
    MONTH_D_CURRENT <- MONTH_D[conditional_selection]
    YEAR_D_CURRENT <- YEAR_D[conditional_selection]
    MONTH_DAY_D_CURRENT <- MONTH_DAY_D[conditional_selection,]

    thresh1 <- .3
    extreme_quantile <- 0.8
    thresh2 <- array(NA,12)
    for (m in 1:12) {
      x <- which(MONTH_D_CURRENT==month_list[m] & PRCP_CURRENT>thresh1)
      thresh2[m] <- quantile(PRCP_CURRENT[x],(extreme_quantile))[[1]]
    }

    #Fit transition probabilities for Markov Chain
    PRCP_LAG0 <- PRCP_CURRENT[2:length(PRCP_CURRENT)]
    PRCP_LAG1 <- PRCP_CURRENT[1:(length(PRCP_CURRENT)-1)]
    MONTH_LAG0 <- MONTH_D_CURRENT[2:length(PRCP_CURRENT)]
    MONTH_LAG1 <- MONTH_D_CURRENT[1:(length(PRCP_CURRENT)-1)]
    YEAR_LAG0 <- YEAR_D_CURRENT[2:length(PRCP_CURRENT)]
    YEAR_LAG1 <- YEAR_D_CURRENT[1:(length(PRCP_CURRENT)-1)]
    for (m in 1:12) {
      x <- which(MONTH_LAG1==month_list[m])
      r <- which(MONTH_SIM==month_list[m] & WATER_YEAR_SIM==(y+START_YEAR_SIM))
      CUR_PRCP0 <- PRCP_LAG0[x]
      CUR_PRCP1 <- PRCP_LAG1[x]
      p00_final[r] <- length(which(PRCP_LAG1[x]<=thresh1 & PRCP_LAG0[x]<=thresh1)) / length(which(PRCP_LAG1[x]<=thresh1))
      p01_final[r] <- length(which(PRCP_LAG1[x]<=thresh1 & PRCP_LAG0[x]>thresh1 & PRCP_LAG0[x]<=thresh2[m])) / length(which(PRCP_LAG1[x]<=thresh1))
      p02_final[r] <- length(which(PRCP_LAG1[x]<=thresh1 & PRCP_LAG0[x]>thresh2[m])) / length(which(PRCP_LAG1[x]<=thresh1))
      p10_final[r] <- length(which(PRCP_LAG1[x]>thresh1 & PRCP_LAG1[x]<=thresh2[m] & PRCP_LAG0[x]<=thresh1)) / length(which(PRCP_LAG1[x]>thresh1 & PRCP_LAG1[x]<=thresh2[m]))
      p11_final[r] <- length(which(PRCP_LAG1[x]>thresh1 & PRCP_LAG1[x]<=thresh2[m] & PRCP_LAG0[x]>thresh1 & PRCP_LAG0[x]<=thresh2[m])) / length(which(PRCP_LAG1[x]>thresh1 & PRCP_LAG1[x]<=thresh2[m]))
      p12_final[r] <- length(which(PRCP_LAG1[x]>thresh1 & PRCP_LAG1[x]<=thresh2[m] & PRCP_LAG0[x]>thresh2[m])) / length(which(PRCP_LAG1[x]>thresh1 & PRCP_LAG1[x]<=thresh2[m]))
      p20_final[r] <- length(which(PRCP_LAG1[x]>thresh2[m] & PRCP_LAG0[x]<=thresh1)) / length(which(PRCP_LAG1[x]>thresh2[m]))
      p21_final[r] <- length(which(PRCP_LAG1[x]>thresh2[m] & PRCP_LAG0[x]>thresh1 & PRCP_LAG0[x]<=thresh2[m])) / length(which(PRCP_LAG1[x]>thresh2[m]))
      p22_final[r] <- length(which(PRCP_LAG1[x]>thresh2[m] & PRCP_LAG0[x]>thresh2[m])) / length(which(PRCP_LAG1[x]>thresh2[m]))

      PI_old <- GET_PI(p00_final[r][1],p01_final[r][1],p02_final[r][1],p10_final[r][1],p11_final[r][1],p12_final[r][1],p20_final[r][1],p21_final[r][1],p22_final[r][1])

      #adjustments for dry spells
      p01_new <- (p01_final[r] + p02_final[r])/dry_spell_change[m] - p02_final[r]
      p00_new <- p00_final[r] + (p01_final[r] - p01_new)
      p01_final[r] <- p01_new
      p00_final[r] <- p00_new

      #adjustments for wet spells
      p10_new <- (p10_final[r] + p12_final[r])/wet_spell_change[m] - p12_final[r]
      p11_new <- p11_final[r] + (p10_final[r] - p10_new)
      p10_final[r] <- p10_new
      p11_final[r] <- p11_new

      PI_new <- GET_PI(p00_final[r][1],p01_final[r][1],p02_final[r][1],p10_final[r][1],p11_final[r][1],p12_final[r][1],p20_final[r][1],p21_final[r][1],p22_final[r][1])
      Ratio_PI_old_PI_new[y,m] <- sum(PI_old[2:3])/sum(PI_new[2:3])
      if (is.na(Ratio_PI_old_PI_new[y,m])) {Ratio_PI_old_PI_new[y,m] <- 1}
    }


    for (j in 1:365) {
      count <- count + 1
      if (count <=SIM_LENGTH) {

        #MARKOV CHAIN
        rn <- rn_all[(count-1)]
        if (OCCURENCES[(count-1)]==0) {
          pp1 <- p00_final[(count-1)]
          pp2 <- p00_final[(count-1)] + p01_final[(count-1)]
        }
        if (OCCURENCES[(count-1)]==1) {
          pp1 <- p10_final[(count-1)]
          pp2 <- p10_final[(count-1)] + p11_final[(count-1)]
        }
        if (OCCURENCES[(count-1)]==2) {
          pp1 <- p20_final[(count-1)]
          pp2 <- p20_final[(count-1)] + p21_final[(count-1)]
        }
        if(rn < pp1) {
          OCCURENCES[count] <- 0
        } else if (rn >= pp1 & rn < pp2) {
          OCCURENCES[count] <- 1
        } else {
          OCCURENCES[count] <- 2
        }

        #KNN RESAMPLING
        m <- MONTH_SIM[(count-1)]
        mmm <- which(month_list==m)
        d <- DAY_SIM[(count-1)]
        cur_OCCERENCE <- OCCURENCES[(count-1)]
        next_OCCURENCE <- OCCURENCES[(count)]

        cur_day <- which(MONTH_DAY_D_CURRENT[,1]==m & MONTH_DAY_D_CURRENT[,2]==d)
        cur_day <- c((cur_day-3),(cur_day-2),(cur_day-1),cur_day,(cur_day+1),(cur_day+2),(cur_day+3))
        cur_day <- subset(cur_day,cur_day > 0)
        if (cur_OCCERENCE==0 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh1)}
        if (cur_OCCERENCE==0 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=thresh1 & PRCP_CURRENT[(cur_day+1)]>thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh2[mmm])}
        if (cur_OCCERENCE==0 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=thresh1 & PRCP_CURRENT[(cur_day+1)]>thresh2[mmm])}
        if (cur_OCCERENCE==1 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh1 & PRCP_CURRENT[cur_day]<=thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]<=thresh1)}
        if (cur_OCCERENCE==1 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh1 & PRCP_CURRENT[cur_day]<=thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh2[mmm])}
        if (cur_OCCERENCE==1 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh1 & PRCP_CURRENT[cur_day]<=thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh2[mmm])}
        if (cur_OCCERENCE==2 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]<=thresh1)}
        if (cur_OCCERENCE==2 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh2[mmm])}
        if (cur_OCCERENCE==2 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh2[mmm])}

        if (length(cur_day_cur_state)==0) {
          cur_day <- which(MONTH_DAY_D_CURRENT[,1]==m & MONTH_DAY_D_CURRENT[,2]==d)
          cur_day_final <- array(NA,length(cur_day)*61)
          cur_day_window <- seq(-30,30)
          for (cc in 1:61) {
            cur_day_final[(1 + length(cur_day)*(cc-1)):(length(cur_day) + length(cur_day)*(cc-1))] <- (cur_day+cur_day_window[cc])
          }
          cur_day <- subset(cur_day_final,cur_day_final > 0)
          if (cur_OCCERENCE==0 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh1)}
          if (cur_OCCERENCE==0 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=thresh1 & PRCP_CURRENT[(cur_day+1)]>thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh2[mmm])}
          if (cur_OCCERENCE==0 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]<=thresh1 & PRCP_CURRENT[(cur_day+1)]>thresh2[mmm])}
          if (cur_OCCERENCE==1 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh1 & PRCP_CURRENT[cur_day]<=thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]<=thresh1)}
          if (cur_OCCERENCE==1 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh1 & PRCP_CURRENT[cur_day]<=thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh2[mmm])}
          if (cur_OCCERENCE==1 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh1 & PRCP_CURRENT[cur_day]<=thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh2[mmm])}
          if (cur_OCCERENCE==2 & next_OCCURENCE==0) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]<=thresh1)}
          if (cur_OCCERENCE==2 & next_OCCURENCE==1) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh1 & PRCP_CURRENT[(cur_day+1)]<=thresh2[mmm])}
          if (cur_OCCERENCE==2 & next_OCCURENCE==2) {cur_day_cur_state <- which(PRCP_CURRENT[cur_day]>thresh2[mmm] & PRCP_CURRENT[(cur_day+1)]>thresh2[mmm])}
        }


        possible_days <- cur_day[cur_day_cur_state]
        PRCP_TODAY <- PRCP_CURRENT[possible_days]
        TEMP_TODAY <- TEMP_CURRENT[possible_days]
        PRCP_TOMORROW <- PRCP_CURRENT[possible_days+1]
        TEMP_TOMORROW <- TEMP_CURRENT[possible_days+1]
        TMAX_TOMORROW <- TMAX_CURRENT[possible_days+1]
        TMIN_TOMORROW <- TMIN_CURRENT[possible_days+1]
        DATE_TOMORROW <- DATE_D_CURRENT[possible_days+1]

        cur_sim_PRCP <- SIM_PRCP[(count-1)]
        cur_sim_TEMP <- SIM_TEMP[(count-1)]

        mm <- which(MONTH_D_CURRENT==m)
        mm_p <- which(MONTH_D_CURRENT==m & PRCP_CURRENT>0)
        sd_monthly_TEMP <- sd(TEMP_CURRENT[mm],na.rm=TRUE)
        sd_monthly_PRCP <- sd(PRCP_CURRENT[mm_p],na.rm=TRUE)
        mean_monthly_TEMP <- mean(TEMP_CURRENT[mm],na.rm=TRUE)
        mean_monthly_PRCP <- mean(PRCP_CURRENT[mm_p],na.rm=TRUE)

        k <- round(sqrt(length(possible_days)))
        RESULT <- KNN(cur_sim_PRCP,cur_sim_TEMP,PRCP_TODAY,TEMP_TODAY,PRCP_TOMORROW,TEMP_TOMORROW,TMAX_TOMORROW,TMIN_TOMORROW,DATE_TOMORROW,k,sd_monthly_PRCP,sd_monthly_TEMP,mean_monthly_PRCP,mean_monthly_TEMP,k1,count)

        SIM_PRCP[count] <- RESULT[1]
        SIM_TEMP[count] <- RESULT[2]
        SIM_TMAX[count] <- RESULT[3]
        SIM_TMIN[count] <- RESULT[4]
        SIM_DATE[count] <- DATE_D_CURRENT[which(as.numeric(DATE_D_CURRENT)==RESULT[5])][1]
      }
    }
  }
  Monthly_Ratio_PI_old_PI_new <- apply(Ratio_PI_old_PI_new,FUN=mean,2)
  return(list(SIM_DATE,YEAR_SIM,MONTH_SIM,DAY_SIM,SIM_PRCP,SIM_TEMP,SIM_TMAX,SIM_TMIN,Monthly_Ratio_PI_old_PI_new))
}
