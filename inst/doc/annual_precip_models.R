## ----load data, warning=FALSE, message=FALSE-----------------------------
library(knitr)
knitr::opts_chunk$set(eval=FALSE)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())
library(weatherGen)
data(climate)

# compute daily averages across sites
clim.da <- group_by(climate, DATE) %>%
  summarise(N=n(),
            PRCP=mean(PRCP),
            TMIN=mean(TMIN),
            TMAX=mean(TMAX))

# aggregate by water year
clim.wyr <- clim.da %>%
  mutate(WYEAR=ifelse(month(DATE)>=10, year(DATE)+1, year(DATE))) %>%
  group_by(WYEAR) %>%
  summarise(N=n(),
            PRCP=sum(PRCP),
            TMIN=mean(TMIN),
            TMAX=mean(TMAX))

## ----plot annual data, fig.cap="Annual Timeseries"-----------------------
#  ggplot(clim.wyr, aes(WYEAR, PRCP)) +
#    geom_line() +
#    labs(x="Water Year", y="Annual Precipitation (mm)")

## ----arima model---------------------------------------------------------
#  arima.models <- list(PRECIP=auto.arima(clim.wyr$PRCP, max.p=2, max.q=2,
#                                         max.P=0, max.Q=0, stationary=TRUE))
#  arima.models.no.warm[["PRECIP"]]

## ----plot_arima----------------------------------------------------------
#  plot(forecast(arima.models.no.warm[["PRECIP"]], h=nrow(clim.wyr)))

## ----, eval=FALSE--------------------------------------------------------
#  wave.spectrum <- wavelet_analysis[['spectrum']]
#  wave.signif <- wavelet_analysis[['signif']]
#  wave.period <- wavelet_analysis[['period']]
#  wave.n <- length(wave.spectrum)
#  
#  sim.n.wyr <- 40
#  
#  iter <- 100
#  
#  sim.prcp.wyr <- array(NA, c(sim.n.wyr, iter))
#  sim.prcp.spectrum <- array(NA, c(wave.n, iter))  # POWER_SPECTRUM_PRCP_ARIMA_SIM

## ----eval=FALSE----------------------------------------------------------
#  for (i in 1:iter) {
#    iter.sim.prcp.wyr <- array(NA, c(sim.n.wyr, length(arima.models)))
#    colnames(iter.sim.prcp.wyr) <- names(arima.models)
#    for (m in names(arima.models)) {
#      cur.model <- arima.models[[m]]
#      cur.model.ar <- as.vector(cur.model$coef)[which(names(cur.model$coef) != "intercept" &
#                                                substr(names(cur.model$coef),1,2)!="ma")]
#      cur.model.ma <- as.vector(cur.model$coef)[which(names(cur.model$coef) != "intercept" &
#                                                substr(names(cur.model$coef),1,2)!="ar")]
#      #in auto.arima, the "intercept" is actually the mean (see http://www.stat.pitt.edu/stoffer/tsa2/Rissues.htm)
#      if ('intercept' %in% names(cur.model$coef)) {
#        # if arima model has an intercept, get it
#        intercept <- cur.model$coef['intercept']
#      } else {
#        intercept <- 0
#      }
#      iter.sim.prcp.wyr[, m] <- arima.sim(n = sim.n.wyr,
#                                          list(ar = cur.model.ar, ma = cur.model.ma),
#                                          sd = sqrt(cur.model$sigma2[[1]])) + intercept
#    }
#    if (ncol(iter.sim.prcp.wyr) > 1) {
#      sim.prcp.wyr[,i] <- apply(iter.sim.prcp.wyr, FUN=sum, 1)
#    } else {
#      sim.prcp.wyr[,i] <- iter.sim.prcp.wyr[, 1]
#    }
#  
#  #   if (transform_data) {
#  #     sim.prcp.wyr[,i] <- InverseBoxCox(lambdas_ANNUAL_PRCP,sim.prcp.wyr[,i])
#  #   }
#  
#    i_wavelet_analysis <- wavelet_annual(sim.prcp.wyr[,i],
#                                         sig=0.90,
#                                         noise.type="white",
#                                         plot.flag=FALSE)
#    sim.prcp.spectrum[,i] <- i_wavelet_analysis[['spectrum']][1:wave.n]
#  }
#  
#  sim.prcp.ts.stat <- cbind(MEAN=rowMeans(sim.prcp.wyr),
#                            Q025=apply(sim.prcp.wyr, 1, quantile, 0.025),
#                            Q975=apply(sim.prcp.wyr, 1, quantile, 0.975))
#  
#  sim.prcp.spectrum.stat <- cbind(MEAN=apply(sim.prcp.spectrum, 1, FUN=mean),
#                                  Q025=apply(sim.prcp.spectrum, 1, FUN=quantile, 0.025),
#                                  Q975=apply(sim.prcp.spectrum, 1, FUN=quantile, 0.975))

## ----plot_sim_precip, eval=FALSE, echo=FALSE-----------------------------
#  par(mfrow=c(1,1))
#  plot(sim.prcp.ts.stat[, 'MEAN'], type='n',
#       xlab="Water Year", ylab="Simulated Precip (mm/yr)",
#       ylim=range(sim.prcp.ts.stat))
#  polygon(c(seq(1, nrow(sim.prcp.ts.stat)),
#            rev(seq(1, nrow(sim.prcp.ts.stat)))),
#          c(sim.prcp.ts.stat[, 'Q025'],
#            rev(sim.prcp.ts.stat[, 'Q975'])),col="grey")
#  lines(1:nrow(sim.prcp.ts.stat), sim.prcp.ts.stat[, 'Q025'])
#  lines(1:nrow(sim.prcp.ts.stat), sim.prcp.ts.stat[, 'Q975'])
#  lines(sim.prcp.ts.stat[, 'MEAN'], col='red')

## ----plot_sim_precip_ggplot, eval=FALSE----------------------------------
#  library(tidyr)
#  as.data.frame(sim.prcp.wyr) %>%
#    mutate(YEAR=seq(1:nrow(sim.prcp.wyr))) %>%
#    gather(TRIAL, VALUE, -YEAR) %>%
#    ggplot(aes(YEAR, VALUE)) +
#    geom_line(aes(group=TRIAL), alpha=0.1) +
#    stat_summary(fun.y=mean, geom='line', color='red') +
#    labs(x='Year', y="Simulated Annual Precip (mm/yr)")

## ----plot_sim_spectrum, eval=FALSE---------------------------------------
#  par(mfrow=c(2,2))
#  x1 <- 1
#  x2 <- wave.n
#  ymin <- min(wave.spectrum[x1:x2],
#              sim.prcp.spectrum.stat[x1:x2, 'MEAN'],
#              wave.signif[x1:x2],
#              sim.prcp.spectrum.stat[x1:x2, 'Q025'])
#  ymax <- max(wave.spectrum[x1:x2],
#              sim.prcp.spectrum.stat[x1:x2, 'MEAN'],
#              wave.signif[x1:x2],
#              sim.prcp.spectrum.stat[x1:x2, 'Q975'])
#  plot(wave.period[x1:x2], wave.spectrum[x1:x2],
#       type="n", ylim=c(ymin,ymax), xlab="Period (years)", ylab="",
#       main="", log="y")
#  mtext(side=2, expression(paste("Power (",mm^2,")")), line=2.5)
#  polygon(c(wave.period[x1:x2],
#            rev(wave.period[x1:x2])),
#          c(sim.prcp.spectrum.stat[x1:x2, 'Q025'],
#            rev(sim.prcp.spectrum.stat[x1:x2, 'Q975'])),col="grey")
#  lines(wave.period[x1:x2], wave.spectrum[x1:x2])
#  lines(wave.period[x1:x2], sim.prcp.spectrum.stat[x1:x2, 'MEAN'], lty=2, col="blue")
#  lines(wave.period[x1:x2], wave.signif[x1:x2], col="red", lty=3)
#  
#  boxplot(apply(sim.prcp.wyr, 2, mean), main="Mean")
#  points(mean(clim.wyr$PRCP),col="red",pch=19)
#  points(mean(apply(sim.prcp.wyr, 2, mean)), col="blue", pch=19)
#  
#  boxplot(apply(sim.prcp.wyr, 2, sd), main="Standard Deviation")
#  points(sd(clim.wyr$PRCP), col="red", pch=19)
#  points(mean(apply(sim.prcp.wyr, 2, sd)), col="blue", pch=19)
#  
#  boxplot(apply(sim.prcp.wyr, 2, skewness), main="Skew")
#  points(skewness(clim.wyr$PRCP), col="red", pch=19)
#  points(mean(apply(sim.prcp.wyr, 2, skewness)),col="blue",pch=19)
#  
#  par(mfrow=c(1,1))

## ----, eval=FALSE--------------------------------------------------------
#  x <- clim.wyr$PRCP
#  n.periods <- 1
#  sig.periods <- seq(8,15)
#  n.comp.periods <- 8
#  
#  scale
#  
#  
#  freq.comps <- array(0,c(length(x),n.periods))
#  
#  for (i in 1:n.periods) {
#    cur.periods <- sig.periods[1:n.comp.periods[i]]
#    if (i>1) {
#      cur.periods <- sig.periods[(1 + (i-1)*n.comp.periods[i-1]):(n.comp.periods[i] + (i-1)*n.comp.periods[i-1])]
#    }
#    sj <- scale[cur.periods]
#  
#    #for Morlet Wavelet with freq = 6
#    Cdelta <- .776
#    w0_0 <- pi^(-1/4)
#  
#    if (length(cur.periods)>1) {
#      freq.comps[,i] <- apply(sd(x)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[cur.periods,]/sqrt(sj),FUN=sum,c(2))
#    }
#    if (length(cur.periods)==1) {
#      freq.comps[,i] <- sd(x)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[cur.periods,]/sqrt(sj)
#    }
#  }
#  
#  noise <- x - apply(freq.comps, 1, FUN=sum)
#  
#  par(mfrow=c((2+n.periods),1))
#  plot(x,type="l",main="ORIGINAL DATA",xlab="TIME (YEARS)",ylab="")
#  for (k in 1:n.periods){
#    plot(freq.comps[,1],type="l",main=paste("COMPONENT",k),xlab="TIME (YEARS)",ylab="")
#  }
#  plot(noise,type="l",main="NOISE",xlab="TIME (YEARS)",ylab="")
#  par(mfrow=c(1,1))
#  

## ----, echo=FALSE, eval=FALSE--------------------------------------------
#  
#  #Change in Lag-1 Autocorrelation of Annual Series (multiplicative)
#  AUTOCOR_ANNUAL_CHANGES <- c(1)
#  
#  #Change in Standard Deviation of Annual Series (multiplicative)
#  MAGNITUDE_NOISE_ANNUAL_CHANGES <- c(1)
#  
#  #Specify Changes by Month (multiplicative except temperature)
#  #Wet Spell Length    	J   F   M   A   M   J   J   A   S   O   N   D
#  E_WET_SPELL_CHANGES1 <- c( .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5)
#  E_WET_SPELL_CHANGES2 <- c(  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
#  E_WET_SPELL_CHANGES3 <- c(  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2)
#  E_WET_SPELL_CHANGES <- rbind(E_WET_SPELL_CHANGES2)
#  
#  #Dry Spell Length		 	J   F   M   A   M   J   J   A   S   O   N   D
#  E_DRY_SPELL_CHANGES1 <- c( .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5)
#  E_DRY_SPELL_CHANGES2 <- c(  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
#  E_DRY_SPELL_CHANGES3 <- c(  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2)
#  E_DRY_SPELL_CHANGES <- rbind(E_DRY_SPELL_CHANGES2)
#  
#  #Mean Prcp		          J   F   M   A   M   J   J   A   S   O   N   D
#  PRCP_MEAN_CHANGES1 <- c( .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9)
#  PRCP_MEAN_CHANGES2 <- c(  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
#  PRCP_MEAN_CHANGES3 <- c(  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1)
#  PRCP_MEAN_CHANGES <- rbind(PRCP_MEAN_CHANGES2)
#  
#  #CV Prcp		        J   F   M   A   M   J   J   A   S   O   N   D
#  PRCP_CV_CHANGES1 <- c( .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9)
#  PRCP_CV_CHANGES2 <- c(  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
#  PRCP_CV_CHANGES3 <- c(  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1)
#  PRCP_CV_CHANGES <- rbind(PRCP_CV_CHANGES2)
#  
#  #Mean Temperature    J   F   M   A   M   J   J   A   S   O   N   D
#  TEMP_CHANGES1 <- c( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
#  TEMP_CHANGES2 <- c(  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
#  TEMP_CHANGES3 <- c(  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
#  TEMP_CHANGES <- rbind(TEMP_CHANGES2)
#  
#  
#  NUM_AUTOCOR_ANNUAL_CHANGES <- length(AUTOCOR_ANNUAL_CHANGES)
#  NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES <- length(MAGNITUDE_NOISE_ANNUAL_CHANGES)
#  NUM_E_WET_SPELL_CHANGES <- dim(E_WET_SPELL_CHANGES)[1]
#  NUM_E_DRY_SPELL_CHANGES <- dim(E_DRY_SPELL_CHANGES)[1]
#  NUM_PRCP_MEAN_CHANGES <- dim(PRCP_MEAN_CHANGES)[1]
#  NUM_PRCP_CV_CHANGES <- dim(PRCP_CV_CHANGES)[1]
#  NUM_TEMP_CHANGES <- dim(TEMP_CHANGES)[1]
#  
#  NEW_AR_MODEL <- list()
#  count <- 0
#  for (k2 in 1:NUM_AUTOCOR_ANNUAL_CHANGES) {
#    for (k3 in 1:NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES) {
#  		count <- count + 1
#  		cur_dir <- main_dir
#  		Lag1_Change <- AUTOCOR_ANNUAL_CHANGES[k2]
#  		Stdev_Change <- MAGNITUDE_NOISE_ANNUAL_CHANGES[k3]
#  		CUR_MODEL <- AR_MODELS[[1]]
#  		e_var_org <- CUR_MODEL$sigma2
#  		AR_org <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
#  		MA_org <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]
#  		Other_e_var <- list()
#  		Other_AR <- list()
#  		Other_MA <- list()
#  		if (USE_WARM_MODEL) {
#  			for (l in 2:dim(Wavelet_Decomposition)[2]) {
#  				CUR_MODEL <- AR_MODELS[[l]]
#  				Other_e_var[[l-1]] <- CUR_MODEL$sigma2
#  				Other_AR[[l-1]] <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
#  				Other_MA[[l-1]] <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]
#  			}
#  		}
#  		new_parameters <- Adjust_WARM_Model(Stdev_Change,Lag1_Change,e_var_org,AR_org,MA_org,Other_e_var,Other_AR,Other_MA,cur_dir)
#  		NEW_AR_MODEL[[count]] <- AR_MODELS[[1]]
#  		NEW_AR_MODEL[[count]]$sigma2 <- new_parameters[1]
#  		non_intercept <- which(names(AR_MODELS[[1]]$coef)!="intercept")
#  		if (length(non_intercept)>0) {
#  			for (cur_par in 1:length(non_intercept)) {NEW_AR_MODEL[[count]]$coef[[cur_par]] <- new_parameters[[cur_par+1]]}
#  		}
#  	}
#  }

## ----, eval=FALSE, echo=FALSE--------------------------------------------
#  
#  setwd(stress_test_dir)
#  unlink(paste(stress_test_dir,"/*",sep=""),recursive=TRUE)
#  Final_Annual_Sim_All <- array(NA,c(num_year_sim,NUM_TRIALS,NUM_AUTOCOR_ANNUAL_CHANGES,NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES))
#  
#  seed_number <- sample(1:10000,size=NUM_TRIALS)
#  
#  for (k1 in 1:NUM_TRIALS) {
#    dir_ending1 <- paste("TRIAL_",k1,sep="")
#  	dir_name1 <- paste(stress_test_dir,dir_ending1,sep="/")
#  	dir.create(dir_name1, showWarnings = TRUE, recursive = FALSE)
#  	setwd(dir_name1)
#  	count <- 0
#  	for (k2 in 1:NUM_AUTOCOR_ANNUAL_CHANGES) {
#  		for (k3 in 1:NUM_MAGNITUDE_NOISE_ANNUAL_CHANGES) {
#  			count <- count + 1
#  			set.seed(seed_number[k1])
#  			Annual_Sim <- array(NA,c(num_year_sim,length(AR_MODELS)))
#  			for (l in 1:length(AR_MODELS)) {
#  				CUR_MODEL <- AR_MODELS[[l]]
#  				if (l==1) {CUR_MODEL <- NEW_AR_MODEL[[count]]}
#  				AR <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ma")]
#  				MA <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)!="intercept" & substr(names(CUR_MODEL$coef),1,2)!="ar")]
#  				INTERCEPT <- 0		#in auto.arima, the "intercept" is actually the mean (see http://www.stat.pitt.edu/stoffer/tsa2/Rissues.htm)
#  				if (length(which(names(CUR_MODEL$coef)=="intercept"))>0) {INTERCEPT <- as.vector(CUR_MODEL$coef)[which(names(CUR_MODEL$coef)=="intercept")]}
#  				Annual_Sim[,l] <- arima.sim(n = num_year_sim, list(ar = AR, ma=MA),sd = sqrt(CUR_MODEL$sigma2[[1]])) + INTERCEPT
#  			}
#  			if (dim(Annual_Sim)[2]>1) {PRCP_FINAL_ANNUAL_SIM <- apply(Annual_Sim,FUN=sum,1)} else {PRCP_FINAL_ANNUAL_SIM <- Annual_Sim[,1]}
#  			Final_Annual_Sim_All[,k1,k2,k3] <- PRCP_FINAL_ANNUAL_SIM
#  			for (k4 in 1:NUM_E_WET_SPELL_CHANGES) {
#  				for (k5 in 1:NUM_E_DRY_SPELL_CHANGES) {
#  					#run weather generator
#  					wet_spell_change <- E_WET_SPELL_CHANGES[k4,]
#  					dry_spell_change <- E_DRY_SPELL_CHANGES[k5,]
#  					WGEN_output <- DAILY_WEATHER_GENERATOR(k1,num_year_sim,PRCP_FINAL_ANNUAL_SIM,ANNUAL_PRCP,WATER_YEAR_A,WATER_YEAR_D,PRCP,TEMP,TMAX,TMIN,DATE_D,MONTH_D,YEAR_D,MONTH_DAY_D,wet_spell_change,dry_spell_change,month_list,water_year_start,water_year_end)
#  
#  					SIM_DATE <- WGEN_output[[1]]
#  					YEAR_SIM <- WGEN_output[[2]]
#  					MONTH_SIM <- WGEN_output[[3]]
#  					DAY_SIM <- WGEN_output[[4]]
#  					if (Markov_Adjust) {Monthly_Ratio_PI_old_PI_new <- WGEN_output[[9]]} else {Monthly_Ratio_PI_old_PI_new <- rep(1,12)}
#  					tt <- match(SIM_DATE,DATE_D)
#  					ANNUAL_BASIN_PRCP_SIM_TEMPORARY <- (aggregate(apply(SITE_PRCP[tt,],1,mean),FUN=sum,by=list(WATER_YEAR_SIM))[,2])
#  					RATIOS <- PRCP_FINAL_ANNUAL_SIM/ANNUAL_BASIN_PRCP_SIM_TEMPORARY
#  					RATIOS[which(RATIOS>(1+max_adjust))] <- (1+max_adjust)
#  					RATIOS[which(RATIOS<(1-max_adjust))] <- (1-max_adjust)
#  
#  					for (k6 in 1:NUM_PRCP_MEAN_CHANGES) {
#  						for (k7 in 1:NUM_PRCP_CV_CHANGES) {
#  							for (k8 in 1:NUM_TEMP_CHANGES) {
#  								dir_ending2 <- paste("TRIAL",k1,"A.P.COR",k2,"A.P.MAG",k3,"D.P.WETSPELL",k4,"D.P.DRYSPELL",k5,"D.P.MEAN",k6,"D.P.CV",k7,"D.T.MEAN",k8,sep="_")
#  								dir_name2 <- paste(dir_name1,dir_ending2,sep="/")
#  								dir.create(dir_name2, showWarnings = TRUE, recursive = FALSE)
#  								setwd(dir_name2)
#  
#  								for (ss in 1:num_site) {
#  									TMAX_SIM_FINAL <- SITE_TMAX[tt,ss]
#  									TMIN_SIM_FINAL <- SITE_TMIN[tt,ss]
#  									if (wind_included) {WIND_SIM_FINAL <- SITE_WIND[tt,ss]}
#  									#annual precipitation adjustment
#  									CUR_SITE_PRCP <- SITE_PRCP[tt,ss]
#  									if (annual_adjust) {
#  										for (y in 1:num_year_sim) {
#  											CUR_SITE_PRCP[(1+365*(y-1)):(365+365*(y-1))] <- CUR_SITE_PRCP[(1+365*(y-1)):(365+365*(y-1))]*(matrix(array(1,365),ncol=1)%*%matrix(RATIOS,nrow=1))[,y]
#  										}
#  									}
#  									PRCP_SIM_FINAL <- CUR_SITE_PRCP
#  
#  									#monthly adjustments
#  									for (m in 1:12) {
#  										aa <- which(MONTH_SIM==m)
#  										TMAX_SIM_FINAL[aa] <- TMAX_SIM_FINAL[aa] + TEMP_CHANGES[k8,m]
#  										TMIN_SIM_FINAL[aa] <- TMIN_SIM_FINAL[aa] + TEMP_CHANGES[k8,m]
#  										a <- which(MONTH_SIM==m & PRCP_SIM_FINAL>0)
#  										w <- which(CUR_SITE_PRCP>0 & MONTH_SIM==m)
#  										if (length(a)>.05*length(which(MONTH_SIM==m))) {
#  											probs <- (1:length(PRCP_SIM_FINAL[a]))/(length(PRCP_SIM_FINAL[a])+1)
#  											dist1 <- fitdistr(PRCP_SIM_FINAL[a],"gamma")
#  											shp1 <- dist1$estimate[[1]]
#  											rate1 <- dist1$estimate[[2]]
#  											shp2 <- shp1/(PRCP_CV_CHANGES[k7,m])^2
#  											rate2 <- rate1/((Monthly_Ratio_PI_old_PI_new[m]*PRCP_MEAN_CHANGES[k6,m])*PRCP_CV_CHANGES[k7,m]^2)
#  											prcp_distr <- qgamma(probs,shape=shp1, rate = rate1)
#  											prcp_distr2 <- qgamma(probs,shape=shp2, rate = rate2)
#  											ratios <- prcp_distr2/prcp_distr
#  
#  											for (j in 1:length(w)) {
#  												multiplier <- approx(prcp_distr,ratios,CUR_SITE_PRCP[w[j]])$y
#  												if (is.na(multiplier)==TRUE & CUR_SITE_PRCP[w[j]]>max(prcp_distr)) {multiplier <- ratios[length(ratios)]}
#  												if (is.na(multiplier)==TRUE & CUR_SITE_PRCP[w[j]]<min(prcp_distr)) {multiplier <- ratios[1]}
#  												PRCP_SIM_FINAL[w[j]] <- multiplier*CUR_SITE_PRCP[w[j]]
#  											}
#  										} else {
#  											for (j in 1:length(w)) {
#  												PRCP_SIM_FINAL[w[j]] <- CUR_SITE_PRCP[w[j]]
#  											}
#  										}
#  									}
#  									if (wind_included) {
#  										FINAL_STOCHASTIC_DAILY <- cbind(YEAR_SIM,MONTH_SIM,DAY_SIM,PRCP_SIM_FINAL,TMAX_SIM_FINAL,TMIN_SIM_FINAL,WIND_SIM_FINAL)
#  										colnames(FINAL_STOCHASTIC_DAILY) <- c("YEAR","MONTH","DAY","PRCP","TMAX","TMIN","WIND")
#  									} else {
#  										FINAL_STOCHASTIC_DAILY <- cbind(YEAR_SIM,MONTH_SIM,DAY_SIM,PRCP_SIM_FINAL,TMAX_SIM_FINAL,TMIN_SIM_FINAL)
#  										colnames(FINAL_STOCHASTIC_DAILY) <- c("YEAR","MONTH","DAY","PRCP","TMAX","TMIN")
#  									}
#  									write.table(FINAL_STOCHASTIC_DAILY,filenames[ss],row.names=FALSE,col.names=FALSE,sep="\t")
#  								}
#  							}
#  						}
#  					}
#  				}
#  			}
#  		}
#  	}
#  }
#  

## ----session-------------------------------------------------------------
#  sessionInfo()

