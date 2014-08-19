## ----load data, warning=FALSE, message=FALSE-----------------------------
library(forecast)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())
library(weatherGen)
data(climate)

# set RNG seed for reproducibility
set.seed(1234)

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
ggplot(clim.wyr, aes(WYEAR, PRCP)) +
  geom_line() +
  labs(x="Water Year", y="Annual Precipitation (mm)")

## ----scott functions, echo=FALSE-----------------------------------------
WAVELET_ANALYSIS <- function(CLIMATE_VARIABLE,siglvl,background_noise,plot_flag) {
  if (missing(plot_flag)) {plot_flag <- TRUE}
	if (missing(siglvl)) {siglvl <- 0.90}
	if (missing(background_noise)) {background_noise <- "white"}
	
	#Define Wavelet function used that returns transform
	waveletf <- function(k,s) {
		nn <- length(k)
		k0 <- 6    #nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
		z <- array(1,nn)
		z[which(k<=0)] <- 0
		expnt <- -((s*k - k0)^2/2)*z
		norm <- sqrt(s*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
		daughter <- norm*exp(expnt)
		daughter <- daughter*z
		fourier_factor <- (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
		coi <- fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
		dofmin <- 2  	# Degrees of freedom	
		return(daughter)
	}

	#Define Wavelet function used that returns fourier_factor, cone of influence, and degrees of freedom
	waveletf2 <- function(k,s) {
		nn <- length(k)
		k0 <- 6    #nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
		z <- array(1,nn)
		z[which(k<=0)] <- 0
		expnt <- -((s*k - k0)^2/2)*z
		norm <- sqrt(s*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
		daughter <- norm*exp(expnt)
		daughter <- daughter*z
		fourier_factor <- (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
		coi <- fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
		dofmin <- 2  	# Degrees of freedom	
		return(c(fourier_factor,coi,dofmin))
	}

	# Perform Wavelet
	#....construct time series to analyze, pad if necessary
	CURRENT_CLIMATE_VARIABLE_org <- CLIMATE_VARIABLE
	variance1 <- var(CURRENT_CLIMATE_VARIABLE_org)
	n1 <- length(CURRENT_CLIMATE_VARIABLE_org)
	CURRENT_CLIMATE_VARIABLE <- scale(CURRENT_CLIMATE_VARIABLE_org)
	variance2 <- var(CURRENT_CLIMATE_VARIABLE)
	base2 <- floor(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
	CURRENT_CLIMATE_VARIABLE <- c(CURRENT_CLIMATE_VARIABLE,rep(0,(2^(base2+1)-n1)))
	n <- length(CURRENT_CLIMATE_VARIABLE)

	#Determine parameters for Wavelet analysis
	dt <- 1
	dj <- 0.25
	s0 <- 2*dt
	J <- floor((1/dj)*log((n1*dt/s0),base=2))

	#....construct SCALE array & empty PERIOD & WAVE arrays
	scale <- s0*2^((0:J)*dj)
	period <- scale
	wave <- array(as.complex(0),c(J+1,n))  # define the wavelet array

	#....construct wavenumber array used in transform [Eqn(5)]
	k <- c(1:floor(n/2))
	k <- k*((2.*pi)/(n*dt))
	k <- c(0,k,-rev(k[1:floor((n-1)/2)]))

	f <- fft(CURRENT_CLIMATE_VARIABLE,inverse=FALSE)        #fourier transform of standardized precipitation

	# loop through all scales and compute transform
	for (a1 in 1:(J+1)) {
		daughter <- waveletf(k,scale[a1])
		results <- waveletf2(k,scale[a1])
		fourier_factor <- results[1]
		coi <- results[2]
		dofmin <- results[3]
		wave[a1,] <- fft(f*daughter,inverse=TRUE)/n  # wavelet transform[Eqn(4)]
	}

	period <- fourier_factor*scale
	coi <- coi*dt*c((0.00001),1:((n1+1)/2-1),rev((1:(n1/2-1))),(0.00001))  # COI [Sec.3g]
	wave <- wave[,1:n1]  # get rid of padding before returning
	POWER <- abs(wave)^2
	GWS <- variance1*apply(POWER,FUN=mean,c(1)) #Global Wavelet Spectrum

	# Significance Testing
	# get the appropriate parameters [see Table(2)]
	k0 <- 6
	empir <- c(2,0.776,2.32,0.60)
	dofmin <- empir[1]     # Degrees of freedom with no smoothing
	Cdelta <- empir[2]     # reconstruction factor
	gamma_fac <- empir[3]  # time-decorrelation factor
	dj0 <- empir[4]       # scale-decorrelation factor
	if (background_noise=="white") {lag1 <- 0}    #for red noise background, lag1 autocorrelation = 0.72, for white noise background, lag1 autocorrelation = 0
	if (background_noise=="red") {lag1 <- .72}

	freq <- dt / period   # normalized frequency
	fft_theor <- (1-lag1^2) / (1-2*lag1*cos(freq*2*pi)+lag1^2)  # [Eqn(16)]			
	fft_theor <- fft_theor  # include time-series variance
	dof <- dofmin

	#ENTIRE POWER SPECTRUM
	chisquare <- qchisq(siglvl,dof)/dof
	signif <- fft_theor*chisquare   # [Eqn(18)]
	sig95 <- ((signif))%o%(array(1,n1))  # expand signif --> (J+1)x(N) array
	sig95 <- POWER / sig95         # where ratio > 1, power is significant

	#TIME_AVERAGED (GLOBAL WAVELET SPECTRUM)
	dof <- n1 - scale
	if (length(dof) == 1) {dof <- array(0,(J+1))+dof}
	dof[which(dof < 1)] <- 1
	dof <- dofmin*sqrt(1 + (dof*dt/gamma_fac / scale)^2 )   
	tt <- which(dof < dofmin)
	dof[tt] <- dofmin
	chisquare_GWS <- array(NA,(J+1))
	signif_GWS <- array(NA,(J+1))
	for (a1 in 1:(J+1)) {
		chisquare_GWS[a1] <- qchisq(siglvl,dof[a1])/dof[a1]
		signif_GWS[a1] <- fft_theor[a1]*variance1*chisquare_GWS[a1]
	}
  
  # PLOT
	period_lower_limit <- 0
	sig_periods <- which(GWS>signif_GWS & period > period_lower_limit)
	if (plot_flag) {
		par(mfrow=c(1,2),font.axis=2,font.lab=2,cex.axis=1.1,cex.lab=1.1)
		xx <- 1:n1
		yy <- period
		image(x=xx,y=yy,z=t(log(POWER,base=2)),xlab="Time (Years)",ylab="Fourier Period (Years)",ylim=rev(range(yy)),log="y",col=(heat.colors(12)))
		lines(xx,coi,lty=2)
		contour(x=xx,y=yy,z=t(sig95),add=TRUE,levels=c(-99,1),labels="")

		xmin <- (min(GWS,signif_GWS))
		xmax <- (max(GWS,signif_GWS))
		plot(GWS,yy,type="b",xlim=c(xmin,xmax),xlab="Global Wavelet Spectrum",ylab="Fourier Period (Years)",log="y",ylim=rev(range(yy)))
		lines(signif_GWS,yy,col="red",lty=2)
		if (background_noise=="white") {tt <- paste(siglvl*100,"% confidence \nlevel for white-noise\n spectrum",sep="")}    
		if (background_noise=="red") {tt <- paste(siglvl*100,"% confidence \nlevel for red-noise\n spectrum",sep="")}
		legend(.2*max(GWS),.65*max(yy),tt,col="red",lty=2,box.lty=0,box.lwd=0,cex=.8)
		mtext("SIGNIFICANT PERIODS:",line=2)
		e <- paste(sig_periods[1])
		if (length(sig_periods)==0) {
			mtext("NONE!",line=1)
		} else if (length(sig_periods)==1) {
			mtext(e,line=1)
		} else {
			for (h in 2:length(sig_periods)) {
				e <- paste(e,",",paste(sig_periods[h]),sep="")
			}
			mtext(e,line=1)
		}
	}
	
	return(list(GWS,signif_GWS,period))
}

WAVELET_DECOMPOSITION <- function(CLIMATE_VARIABLE,NUM_FINAL_PERIODS,ALL_SIG_PERIODS,NUM_PERIODS_ALL_COMPS,plot_flag) {
  if (missing(plot_flag)) {plot_flag <- TRUE}
  #Define Wavelet function used that returns transform
	waveletf <- function(k,s) {
		nn <- length(k)
		k0 <- 6    #nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
		z <- array(1,nn)
		z[which(k<=0)] <- 0
		expnt <- -((s*k - k0)^2/2)*z
		norm <- sqrt(s*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
		daughter <- norm*exp(expnt)
		daughter <- daughter*z
		fourier_factor <- (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
		coi <- fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
		dofmin <- 2  	# Degrees of freedom	
		return(daughter)
	}

	#Define Wavelet function used that returns fourier_factor, cone of influence, and degrees of freedom
	waveletf2 <- function(k,s) {
		nn <- length(k)
		k0 <- 6    #nondimensional frequency, here taken to be 6 to satisfy the admissibility condition [Farge 1992]
		z <- array(1,nn)
		z[which(k<=0)] <- 0
		expnt <- -((s*k - k0)^2/2)*z
		norm <- sqrt(s*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
		daughter <- norm*exp(expnt)
		daughter <- daughter*z
		fourier_factor <- (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
		coi <- fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
		dofmin <- 2  	# Degrees of freedom	
		return(c(fourier_factor,coi,dofmin))
	}
	# Perform wavelet transform
	#....construct time series to analyze, pad if necessary
	CURRENT_CLIMATE_VARIABLE_org <- CLIMATE_VARIABLE
	variance1 <- var(CURRENT_CLIMATE_VARIABLE_org)
	n1 <- length(CURRENT_CLIMATE_VARIABLE_org)
	CURRENT_CLIMATE_VARIABLE <- scale(CURRENT_CLIMATE_VARIABLE_org)
	variance2 <- var(CURRENT_CLIMATE_VARIABLE)
	base2 <- floor(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
	CURRENT_CLIMATE_VARIABLE <- c(CURRENT_CLIMATE_VARIABLE,rep(0,(2^(base2+1)-n1)))
	n <- length(CURRENT_CLIMATE_VARIABLE)

	#Determine parameters for Wavelet analysis
	dt <- 1
	dj <- 0.25
	s0 <- 2*dt
	J <- floor((1/dj)*log((n1*dt/s0),base=2))

	#....construct SCALE array & empty PERIOD & WAVE arrays
	scale <- s0*2^((0:J)*dj)
	period <- scale
	wave <- array(as.complex(0),c(J+1,n))  # define the wavelet array

	#....construct wavenumber array used in transform [Eqn(5)]
	k <- c(1:floor(n/2))
	k <- k*((2.*pi)/(n*dt))
	k <- c(0,k,-rev(k[1:floor((n-1)/2)]))

	f <- fft(CURRENT_CLIMATE_VARIABLE,inverse=FALSE)        #fourier transform of standardized precipitation

	# loop through all scales and compute transform
	for (a1 in 1:(J+1)) {
		daughter <- waveletf(k,scale[a1])
		results <- waveletf2(k,scale[a1])
		fourier_factor <- results[1]
		coi <- results[2]
		dofmin <- results[3]
		wave[a1,] <- fft(f*daughter,inverse=TRUE)/n  # wavelet transform[Eqn(4)]
	}

	period <- fourier_factor*scale
	coi <- coi*dt*c((0.00001),1:((n1+1)/2-1),rev((1:(n1/2-1))),(0.00001))  # COI [Sec.3g]
	wave <- wave[,1:n1]  # get rid of padding before returning
	POWER <- abs(wave)^2
	GWS <- variance1*apply(POWER,FUN=mean,c(1)) #Global Wavelet Spectrum

	# Signficance Testing
	background_noise <- "white"   #can be red or white 

	siglvl <- 0.95
	# get the appropriate parameters [see Table(2)]
	k0 <- 6
	empir <- c(2,0.776,2.32,0.60)
	dofmin <- empir[1]     # Degrees of freedom with no smoothing
	Cdelta <- empir[2]     # reconstruction factor
	gamma_fac <- empir[3]  # time-decorrelation factor
	dj0 <- empir[4]       # scale-decorrelation factor
	if (background_noise=="white") {lag1 <- 0}    #for red noise background, lag1 autocorrelation = 0.72, for white noise background, lag1 autocorrelation = 0
	if (background_noise=="red") {lag1 <- .72}

	freq <- dt / period   # normalized frequency
	fft_theor <- (1-lag1^2) / (1-2*lag1*cos(freq*2*pi)+lag1^2)  # [Eqn(16)]			
	fft_theor <- fft_theor  # include time-series variance
	dof <- dofmin

	#ENTIRE POWER SPECTRUM
	chisquare <- qchisq(siglvl,dof)/dof
	signif <- fft_theor*chisquare   # [Eqn(18)]
	sig95 <- ((signif))%o%(array(1,n1))  # expand signif --> (J+1)x(N) array
	sig95 <- POWER / sig95         # where ratio > 1, power is significant

	#TIME_AVERAGED (GLOBAL WAVELET SPECTRUM)
	dof <- n1 - scale
	if (length(dof) == 1) {dof <- array(0,(J+1))+dof}
	dof[which(dof < 1)] <- 1
	dof <- dofmin*sqrt(1 + (dof*dt/gamma_fac / scale)^2 )   
	tt <- which(dof < dofmin)
	dof[tt] <- dofmin
	chisquare_GWS <- array(NA,(J+1))
	signif_GWS <- array(NA,(J+1))
	for (a1 in 1:(J+1)) {
		chisquare_GWS[a1] <- qchisq(siglvl,dof[a1])/dof[a1]
		signif_GWS[a1] <- fft_theor[a1]*variance1*chisquare_GWS[a1]
	}

	LOW_FREQUENCY_COMPONENTS <- array(0,c(length(CLIMATE_VARIABLE),NUM_FINAL_PERIODS))
	for (i in 1:NUM_FINAL_PERIODS) {
		CUR_PERIODS <- ALL_SIG_PERIODS[1:NUM_PERIODS_ALL_COMPS[i]]
		if (i>1) {CUR_PERIODS <- ALL_SIG_PERIODS[(1 + (i-1)*NUM_PERIODS_ALL_COMPS[i-1]):(NUM_PERIODS_ALL_COMPS[i] + (i-1)*NUM_PERIODS_ALL_COMPS[i-1])]}
		sj <- scale[CUR_PERIODS]
		#for Morlet Wavelet with freq = 6
		Cdelta <- .776
		w0_0 <- pi^(-1/4)
		if (length(CUR_PERIODS)>1) {LOW_FREQUENCY_COMPONENTS[,i] <- apply(sd(CURRENT_CLIMATE_VARIABLE_org)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[CUR_PERIODS,]/sqrt(sj),FUN=sum,c(2))}
		if (length(CUR_PERIODS)==1) {LOW_FREQUENCY_COMPONENTS[,i] <- sd(CURRENT_CLIMATE_VARIABLE_org)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[CUR_PERIODS,]/sqrt(sj)}
	}

	NOISE <- CURRENT_CLIMATE_VARIABLE_org - apply(LOW_FREQUENCY_COMPONENTS,FUN=sum,c(1)) 
	
	if (plot_flag) {
		par(mfrow=c((2+NUM_FINAL_PERIODS),1),font.axis=2,font.lab=2)
		plot(CLIMATE_VARIABLE,type="l",main="ORIGINAL DATA",xlab="TIME (YEARS)",ylab="")
		for (k in 1:NUM_FINAL_PERIODS){
			plot(LOW_FREQUENCY_COMPONENTS[,1],type="l",main=paste("COMPONENT",k),xlab="TIME (YEARS)",ylab="")
		}
		plot(NOISE,type="l",main="NOISE",xlab="TIME (YEARS)",ylab="")
    par(mfrow=c(1,1))
	}
	return(cbind(NOISE,LOW_FREQUENCY_COMPONENTS))	
}

## ----wavelet analysis scott----------------------------------------------
wave <- WAVELET_ANALYSIS(clim.wyr$PRCP, siglvl=0.90,
                         background_noise="white", plot_flag=TRUE)

## ----wavelet components scott--------------------------------------------
wave.comp <- WAVELET_DECOMPOSITION(CLIMATE_VARIABLE=clim.wyr$PRCP,
                                   NUM_FINAL_PERIODS=1,
                                   ALL_SIG_PERIODS=c(11,12),
                                   NUM_PERIODS_ALL_COMPS=2,
                                   plot_flag=TRUE)
head(wave.comp)

## ----warm arima models---------------------------------------------------
models <- list(NOISE=auto.arima(wave.comp[, 1],
                                max.p=2,max.q=2,max.P=0,max.Q=0,
                                stationary=TRUE),
               COMPONENT=auto.arima(wave.comp[, 2],
                                    max.p=2,max.q=2,max.P=0,max.Q=0,
                                    stationary=TRUE))
models

## ----plot_warm_forecast--------------------------------------------------
par(mfrow=c(2,1))
plot(forecast(models[["NOISE"]], h=nrow(clim.wyr)))
plot(forecast(models[["COMPONENT"]], h=nrow(clim.wyr)))
par(mfrow=c(1,1))

## ----arima model---------------------------------------------------------
ar_models.arima <- list(PRECIP=auto.arima(clim.wyr$PRCP,
                                          max.p=2,max.q=2,max.P=0,max.Q=0,
                                          stationary=TRUE))
ar_models.arima[["PRECIP"]]

## ----plot_arima----------------------------------------------------------
plot(forecast(ar_models.arima[["PRECIP"]], h=nrow(clim.wyr)))

## ----simulate arima------------------------------------------------------
simulate.arima <- function(model, n) {
  # simulate one arima model
  sim <- arima.sim(n=n, 
                   list(ar=coef(model)[grepl('ar', names(coef(model)))],
                        ma=coef(model)[grepl('ma', names(coef(model)))]),
                   sd = sqrt(model$sigma2[[1]]))
  
  # extract intercept
  if ('intercept' %in% names(model$coef)) {
    intercept <- model$coef['intercept']
  } else {
    intercept <- 0
  }

  # add intercept (mean) to current simulation
  sim <- sim + intercept
  
  return(coredata(sim))
}
simulate.arima(model=models[['COMPONENT']], n=40)

## ----simulate arimas-----------------------------------------------------
simulate.arimas <- function(models, n, components=TRUE) {
  # n: number of simulation timesteps
  # models: list of arima models
  # components: if TRUE returns list of $sum and $components, otherwise returns just sum vector
  
  # run simulation on each model and combine into 2-d array
  sim.components <- lapply(models, simulate.arima, n=40) %>% sapply(cbind)
  
  # compute sum of individual component simulations
  if (ncol(sim.components) > 1) {
    sim.sum <- apply(sim.components, FUN=sum, 1)
  } else {
    sim.sum <- sim.components
  }
  
  if (components) {
    return(list(components=sim.components, sum=sim.sum))  
  } else {
    return(sim.sum)
  }
  
}
sim.warm <- simulate.arimas(models=models, n=40)
sim.warm

## ----plot_sim_prcp-------------------------------------------------------
plot(sim.warm$sum, type='l', xlab='Timestep (year)', ylab='Simulated Annual Precip (mm)')

## ----wavelet analysis of sim precip--------------------------------------
wt.sim.warm <- WAVELET_ANALYSIS(sim.warm$sum, siglvl=0.90,
                                background_noise="white", plot_flag=TRUE)

## ----sim arimas monte carlo----------------------------------------------
n.iter <- 100
n.sim <- 40

sim.warm.prcp <- array(NA, c(n.sim, n.iter))
sim.warm.gws <- array(NA, c(length(wt.sim.warm[[1]]), n.iter))
for (i in 1:n.iter) {
  i.sim.warm <- simulate.arimas(models=models, n=40)
  i.wt <- WAVELET_ANALYSIS(i.sim.warm$sum, siglvl=0.90,
                           background_noise="white", plot_flag=FALSE)
  sim.warm.prcp[,i] <- i.sim.warm$sum
  sim.warm.gws[,i] <- i.wt[[1]]
}

sim.warm.prcp.stat <- cbind(MEAN=rowMeans(sim.warm.prcp),
                            Q025=apply(sim.warm.prcp, 1, quantile, 0.025),
                            Q975=apply(sim.warm.prcp, 1, quantile, 0.975))

sim.warm.gws.stat <- cbind(MEAN=rowMeans(sim.warm.gws),
                           Q025=apply(sim.warm.gws, 1, quantile, 0.025),
                           Q975=apply(sim.warm.gws, 1, quantile, 0.975))

## ----plot_sim_precip-----------------------------------------------------
par(mfrow=c(1,1))
plot(sim.warm.prcp.stat[, 'MEAN'], type='n', 
     xlab="Water Year", ylab="Simulated Precip (mm/yr)",
     ylim=range(sim.warm.prcp.stat))
polygon(c(seq(1, nrow(sim.warm.prcp.stat)),
          rev(seq(1, nrow(sim.warm.prcp.stat)))),
        c(sim.warm.prcp.stat[, 'Q025'],
          rev(sim.warm.prcp.stat[, 'Q975'])),col="grey")
lines(sim.warm.prcp.stat[, 'MEAN'], col='red')

## ----plot_sim_precip_ggplot----------------------------------------------
library(tidyr)
as.data.frame(sim.warm.prcp) %>% 
  mutate(YEAR=seq(1:nrow(sim.warm.prcp))) %>%
  gather(TRIAL, VALUE, -YEAR) %>% 
  ggplot(aes(YEAR, VALUE)) +
  geom_line(aes(group=TRIAL), alpha=0.1) +
  stat_summary(fun.y=mean, geom='line', color='red') +
  labs(x='Year', y="Simulated Annual Precip (mm/yr)")

## ----plot_sim_gws, fig.height=8------------------------------------------
par(mfrow=c(2,2))
x1 <- 1
x2 <- length(wave[[3]])
ymin <- min(wave[[1]], wave[[2]],
            sim.warm.gws.stat[x1:x2, 'MEAN'],
            sim.warm.gws.stat[x1:x2, 'Q025'])
ymax <- max(wave[[1]], wave[[2]],
            sim.warm.gws.stat[x1:x2, 'MEAN'],
            sim.warm.gws.stat[x1:x2, 'Q975'])
plot(wave[[3]], wave[[1]],
     type="n", ylim=c(ymin,ymax), xlab="Period (years)", ylab="",
     main="", log="y")
mtext(side=2, expression(paste("Power (",mm^2,")")), line=2.5)
polygon(c(wave[[3]],
          rev(wave[[3]])),
        c(sim.warm.gws.stat[x1:x2, 'Q025'],
          rev(sim.warm.gws.stat[x1:x2, 'Q975'])),
        col="grey")
lines(wave[[3]], wave[[1]])
lines(wave[[3]], sim.warm.gws.stat[x1:x2, 'MEAN'], lty=2, col="blue")
lines(wave[[3]], wave[[2]], col="red", lty=3)

boxplot(colMeans(sim.warm.prcp), main="Mean")
points(mean(clim.wyr$PRCP),col="red",pch=19)
points(mean(colMeans(sim.warm.prcp)), col="blue", pch=19)

boxplot(apply(sim.warm.prcp, 2, sd), main="Standard Deviation")
points(sd(clim.wyr$PRCP), col="red", pch=19)
points(mean(apply(sim.warm.prcp, 2, sd)), col="blue", pch=19)

boxplot(apply(sim.warm.prcp, 2, skewness), main="Skew")
points(skewness(clim.wyr$PRCP), col="red", pch=19)
points(mean(apply(sim.warm.prcp, 2, skewness)),col="blue",pch=19)

## ----scott sce, echo=FALSE-----------------------------------------------
# source: http://ihacres.googlecode.com/svn-history/r89/trunk/R/sce.R
sceDefaults <- function()
    list(ncomplex = 5, ## number of complexes
         cce.iter = NA, ## number of iteration in inner loop (CCE algorithm)
         fnscale = 1, ## function scaling factor (set to -1 for maximisation)
         elitism = 1, ## controls amount of weighting in sampling towards the better parameter sets
         initsample = "latin", ## sampling scheme for initial values -- "latin" or "random"
         reltol = 1e-5, ## convergence threshold: relative improvement factor required in an SCE iteration
         tolsteps = 7, ## number of iterations within reltol to confirm convergence
         maxit = 10000, ## maximum number of iterations
         maxeval = Inf, ## maximum number of function evaluations
         maxtime = Inf, ## maximum duration of optimization in seconds
         returnpop = FALSE, ## whether to return populations from all iterations
         trace = 0, ## level of user feedback
         REPORT = 1) ## number of iterations between reports when trace >= 1


SCEoptim <- function(FUN, par, ..., lower = -Inf, upper = Inf, control = list()) {
    FUN <- match.fun(FUN)
    stopifnot(is.numeric(par))
    stopifnot(length(par) > 0)
    stopifnot(is.numeric(lower))
    stopifnot(is.numeric(upper))
    ## allow `lower` or `upper` to apply to all parameters
    if (length(lower) == 1)
        lower <- rep(lower, length = length(par))
    if (length(upper) == 1)
        upper <- rep(upper, length = length(par))
    stopifnot(length(lower) == length(par))
    stopifnot(length(upper) == length(par))

    ## determine number of variables to be optimized
    NDIM <- length(par)

    ## update default options with supplied options
    stopifnot(is.list(control))
    control <- modifyList(sceDefaults(), control)
    isValid <- names(control) %in% names(sceDefaults())
    if (any(!isValid))
        stop("unrecognised options: ",
             toString(names(control)[!isValid]))

    returnpop <- control$returnpop
    trace <- control$trace
    nCOMPLEXES <- control$ncomplex
    CCEITER <- control$cce.iter
    MAXIT <- control$maxit
    MAXEVAL <- control$maxeval

    ## recommended number of CCE steps in Duan et al 1994:
    if (is.na(CCEITER))
        CCEITER <- 2 * NDIM + 1

    if (is.finite(MAXEVAL)) {
        ## upper bound on number of iterations to reach MAXEVAL
        MAXIT <- min(MAXIT, ceiling(MAXEVAL / (nCOMPLEXES * CCEITER)))
    }

    ## define number of points in each complex
    nPOINTS_COMPLEX <- 2 * NDIM + 1

    ## define number of points in each simplex
    nPOINTS_SIMPLEX <- NDIM+1

    ## define total number of points
    nPOINTS <- nCOMPLEXES * nPOINTS_COMPLEX

    ## initialize counters
    funevals <- 0


    costFunction <- function(FUN, par, ...)
    {
        ## check lower and upper bounds
        i <- which(par < lower)
        if (any(i)) {
            i <- i[1]
            return( 1e12 + (lower[i] - par[i]) * 1e6 )
        }
        i <- which(par > upper)
        if (any(i)) {
            i <- i[1]
            return( 1e12 + (par[i] - upper[i]) * 1e6 )
        }
        funevals <<- funevals + 1
        result <- FUN(par, ...) * control$fnscale
        if (is.na(result))
            result <- 1e12
        result
    }

    simplexStep <- function(P, FAC)
    {
        ## Extrapolates by a factor FAC through the face of the simplex across from
        ## the highest (i.e. worst) point.
        worst <- nPOINTS_SIMPLEX
        centr <- apply(P[-worst,,drop=FALSE], 2, mean)
        newpar <- centr*(1-FAC) + P[worst,]*FAC
        newpar
    }


    ## initialize population matrix
    POPULATION <- matrix(as.numeric(NA), nrow = nPOINTS, ncol = NDIM)
    if (!is.null(names(par)))
        colnames(POPULATION) <- names(par)
    POP.FITNESS <- numeric(length = nPOINTS)
    POPULATION[1,] <- par

    ## generate initial parameter values by random uniform sampling
    finitelower <- ifelse(is.infinite(lower), -(abs(par)+2)*5, lower)
    finiteupper <- ifelse(is.infinite(upper), +(abs(par)+2)*5, upper)
    if (control$initsample == "latin") {
        for (i in 1:NDIM) {
            tmp <- seq(finitelower[i], finiteupper[i], length = nPOINTS-1)
            tmp <- jitter(tmp, factor = 2)
            tmp <- pmax(finitelower[i], pmin(finiteupper[i], tmp))
            POPULATION[-1,i] <- sample(tmp)
        }
    } else {
        for (i in 1:NDIM)
            POPULATION[-1,i] <- runif(nPOINTS-1, finitelower[i], finiteupper[i])
    }

    ## only store all iterations if requested -- could be big!
    if (!is.finite(MAXIT)) {
        MAXIT <- 10000
        warning("setting maximum iterations to 10000")
    }
    if (returnpop) {
        POP.ALL <- array(as.numeric(NA), dim = c(nPOINTS, NDIM, MAXIT))
        if (!is.null(names(par)))
            dimnames(POP.ALL)[[2]] <- names(par)
    }
    POP.FIT.ALL <- matrix(as.numeric(NA), ncol = nPOINTS, nrow = MAXIT)
    BESTMEM.ALL <- matrix(as.numeric(NA), ncol = NDIM, nrow = MAXIT)
    if (!is.null(names(par)))
        colnames(BESTMEM.ALL) <- names(par)

    ## the output object
    obj <- list()
    class(obj) <- c("SCEoptim", class(obj))
    obj$call <- match.call()
    obj$control <- control

    EXITFLAG <- NA
    EXITMSG <- NULL

    ## initialize timer
    tic <- as.numeric(Sys.time())
    toc <- 0

    ## calculate cost for each point in initial population
    for (i in 1:nPOINTS)
        POP.FITNESS[i] <- costFunction(FUN, POPULATION[i,], ...)

    ## sort the population in order of increasing function values
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx,,drop=FALSE]

    ## store one previous iteration only
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS

    if (returnpop) {
        POP.ALL[,,1] <- POPULATION
    }
    POP.FIT.ALL[1,] <- POP.FITNESS
    BESTMEM.ALL[1,] <- POPULATION[1,]

    ## store best solution from last two iterations
    prevBestVals <- rep(Inf, control$tolsteps)
    prevBestVals[1] <- POP.FITNESS[1]

    ## for each iteration...
    i <- 0
    while (i < MAXIT) {

        i <- i + 1

        ## The population matrix POPULATION will now be rearranged into complexes.

        ## For each complex ...
        for (j in 1:nCOMPLEXES) {

            ## construct j-th complex from POPULATION

            k1 <- 1:nPOINTS_COMPLEX
            k2 <- (k1-1) * nCOMPLEXES + j

            COMPLEX <- POP.PREV[k2,,drop=FALSE]
            COMPLEX_FITNESS <- POP.FIT.PREV[k2]

            ## Each complex evolves a number of steps according to the competitive
            ## complex evolution (CCE) algorithm as described in Duan et al. (1992).
            ## Therefore, a number of 'parents' are selected from each complex which
            ## form a simplex. The selection of the parents is done so that the better
            ## points in the complex have a higher probability to be selected as a
            ## parent. The paper of Duan et al. (1992) describes how a trapezoidal
            ## probability distribution can be used for this purpose.

            for (k in 1:CCEITER) {

                ## select simplex by sampling the complex

                ## sample points with "trapezoidal" i.e. linear probability
                weights <- rev(ppoints(nPOINTS_COMPLEX))
                ## 'elitism' parameter can give more weight to the better results:
                weights <- weights ^ control$elitism
                LOCATION <- sample(seq(1,nPOINTS_COMPLEX), size = nPOINTS_SIMPLEX,
                                   prob = weights)

                LOCATION <- sort(LOCATION)

                ## construct the simplex
                SIMPLEX <- COMPLEX[LOCATION,,drop=FALSE]
                SIMPLEX_FITNESS <- COMPLEX_FITNESS[LOCATION]

                worst <- nPOINTS_SIMPLEX

                ## generate new point for simplex

                ## first extrapolate by a factor -1 through the face of the simplex
                ## across from the high point,i.e.,reflect the simplex from the high point
                parRef <- simplexStep(SIMPLEX, FAC = -1)
                fitRef <- costFunction(FUN, parRef, ...)

                ## check the result
                if (fitRef <= SIMPLEX_FITNESS[1]) {
                    ## gives a result better than the best point,so try an additional
                    ## extrapolation by a factor 2
                    parRefEx <- simplexStep(SIMPLEX, FAC = -2)
                    fitRefEx <- costFunction(FUN, parRefEx, ...)
                    if (fitRefEx < fitRef) {
                        SIMPLEX[worst,] <- parRefEx
                        SIMPLEX_FITNESS[worst] <- fitRefEx
                        ALGOSTEP <- 'reflection and expansion'
                    } else {
                        SIMPLEX[worst,] <- parRef
                        SIMPLEX_FITNESS[worst] <- fitRef
                        ALGOSTEP <- 'reflection'
                    }
                } else if (fitRef >= SIMPLEX_FITNESS[worst-1]) {
                    ## the reflected point is worse than the second-highest, so look
                    ## for an intermediate lower point, i.e., do a one-dimensional
                    ## contraction
                    parCon <- simplexStep(SIMPLEX, FAC = -0.5)
                    fitCon <- costFunction(FUN, parCon, ...)
                    if (fitCon < SIMPLEX_FITNESS[worst]) {
                        SIMPLEX[worst,] <- parCon
                        SIMPLEX_FITNESS[worst] <- fitCon
                        ALGOSTEP <- 'one dimensional contraction'
                    } else {
                        ## can't seem to get rid of that high point, so better contract
                        ## around the lowest (best) point
                        SIMPLEX <- (SIMPLEX + rep(SIMPLEX[1,], each=nPOINTS_SIMPLEX)) / 2
                        for (k in 2:NDIM)
                            SIMPLEX_FITNESS[k] <- costFunction(FUN, SIMPLEX[k,], ...)
                        ALGOSTEP <- 'multiple contraction'
                    }
                } else {
                    ## if better than second-highest point, use this point
                    SIMPLEX[worst,] <- parRef
                    SIMPLEX_FITNESS[worst] <- fitRef
                    ALGOSTEP <- 'reflection'
                }

                if (trace >= 3) {
                    message(ALGOSTEP)
                }

                ## replace the simplex into the complex
                COMPLEX[LOCATION,] <- SIMPLEX
                COMPLEX_FITNESS[LOCATION] <- SIMPLEX_FITNESS

                ## sort the complex
                idx <- order(COMPLEX_FITNESS)
                COMPLEX_FITNESS <- COMPLEX_FITNESS[idx]
                COMPLEX <- COMPLEX[idx,,drop=FALSE]
            }

            ## replace the complex back into the population
            POPULATION[k2,] <- COMPLEX
            POP.FITNESS[k2] <- COMPLEX_FITNESS
        }

        ## At this point, the population was divided in several complexes, each of which
        ## underwent a number of iteration of the simplex (Metropolis) algorithm. Now,
        ## the points in the population are sorted, the termination criteria are checked
        ## and output is given on the screen if requested.

        ## sort the population
        idx <- order(POP.FITNESS)
        POP.FITNESS <- POP.FITNESS[idx]
        POPULATION <- POPULATION[idx,,drop=FALSE]
        if (returnpop) {
            POP.ALL[,,i] <- POPULATION
        }
        POP.FIT.ALL[i,] <- POP.FITNESS
        BESTMEM.ALL[i,] <- POPULATION[1,]

        curBest <- POP.FITNESS[1]

        ## end the optimization if one of the stopping criteria is met

        prevBestVals <- c(curBest, head(prevBestVals, -1))
        reltol <- control$reltol
        if (all(abs(diff(prevBestVals)) <= reltol * (abs(curBest)+reltol))) {
            EXITMSG <- 'Change in solution over [tolsteps] less than specified tolerance (reltol).'
            EXITFLAG <- 0
        }

        ## give user feedback on screen if requested
        if (trace >= 1) {
            if (i == 1) {
                message(' Nr Iter Nr Fun Eval Current best function Current worst function')
            }
            if ((i %% control$REPORT == 0) || (!is.na(EXITFLAG)))
            {
                message(sprintf(' %5.0f %5.0f %12.6g %12.6g',
                        i, funevals, min(POP.FITNESS), max(POP.FITNESS)))
                if (trace >= 2)
                    message("parameters: ", toString(signif(POPULATION[1,], 3)))
            }
        }

        if (!is.na(EXITFLAG))
            break

        if ((i >= control$maxit) || (funevals >= control$maxeval)) {
            EXITMSG <- 'Maximum number of function evaluations or iterations reached.'
            EXITFLAG <- 1
            break
        }

        toc <- as.numeric(Sys.time()) - tic
        if (toc > control$maxtime) {
            EXITMSG <- 'Exceeded maximum time.'
            EXITFLAG <- 2
            break
        }

        ## go to next iteration
        POP.PREV <- POPULATION
        POP.FIT.PREV <- POP.FITNESS
    }
    if (trace >= 1)
        message(EXITMSG)

    ## return solution
    obj$par <- POPULATION[1,]
    obj$value <- POP.FITNESS[1]
    obj$convergence <- EXITFLAG
    obj$message <- EXITMSG

    ## store number of function evaluations
    obj$counts <- funevals
    ## store number of iterations
    obj$iterations <- i
    ## store the amount of time taken
    obj$time <- toc

    if (returnpop) {
        ## store information on the population at each iteration
        obj$POP.ALL <- POP.ALL[,,1:i]
        dimnames(obj$POP.ALL)[[3]] <- paste("iteration", 1:i)
    }
    obj$POP.FIT.ALL <- POP.FIT.ALL[1:i,]
    obj$BESTMEM.ALL <- BESTMEM.ALL[1:i,]

    obj
}


## ----scot arma modify, echo=FALSE----------------------------------------
VAR_COV_ARMA <- function(parameters,num_AR,num_MA) {
  e_var <- parameters[1]
	AR <-NULL
	MA <- NULL
	if (num_AR>0) AR <- parameters[2:(2+num_AR-1)]
	if (num_MA>0) MA <- parameters[(num_AR+2):(1+num_AR+num_MA)]

	#test for nonstationarity
	NonStationary <- FALSE
	if (length(AR)==1) {
		if(abs(AR[1])>=1) {NonStationary<-TRUE}
	} else if (length(AR==2)) {
		a1 <- AR[1]
		a2 <- AR[2]
		d <- a1^2 + 4*a2
		if(d>0) {
			test1 <- a1 + a2
			test2 <- a2-a1
			if (test1>=1 | test2>=1) {
				NonStationary<-TRUE
			}
		} else if (d==0) {
			test1 <- abs(a1)
			if (test1>=2) {
				NonStationary<-TRUE
			}		
		} else {
			test1 <- -a2
			if (test1>=1) {
				
				NonStationary<-TRUE
			}				
		}
	}
	
	a1 <- NULL
	a2 <- NULL
	b1 <- NULL
	b2 <- NULL
	#AR(0)
	if (length(AR)==0 & length(MA)==0 & NonStationary==FALSE) {
		VAR <- e_var
		LAG1 <- 0

	#AR(1)
	} else if (length(AR)==1 & length(MA)==0 & NonStationary==FALSE) {
		a1 <- AR[1]
		VAR <- e_var/(1-a1^2)
		LAG1 <- a1

	#AR(2)
	} else if (length(AR)==2 & length(MA)==0 & NonStationary==FALSE) {
		a1 <- AR[1]
		a2 <- AR[2]
		VAR <- ((1-a2)/(1+a2))*e_var/(a1+a2-1)/(a2-a1-1)
		LAG1 <- a1/(1-a2)

	#MA(1)
	} else if (length(AR)==0 & length(MA)==1 & NonStationary==FALSE) {
		b1 <- MA[1]
		VAR <- (1+b1^2)*e_var
		LAG1 <- b1/(1+b1^2)

	#MA(2)
	} else if (length(AR)==0 & length(MA)==2 & NonStationary==FALSE) {
		b1 <- MA[1]
		b2 <- MA[2]
		VAR <- e_var*(1+b1^2+b2^2)
		LAG1 <- b1*(1+b2)/(1+b1^2+b2^2)

	#ARMA(1,1)
	} else if (length(AR)==1 & length(MA)==1 & NonStationary==FALSE) {
		a1 <- AR[1]
		b1 <- MA[1]
		VAR <- (1+b1^2+2*a1*b1)/(1-a1^2)*e_var
		LAG1 <- (1+a1*b1)*(a1+b1)/(1+b1^2+2*a1*b1)
		
	#ARMA(2,1)
	} else if (length(AR)==2 & length(MA)==1 & NonStationary==FALSE) {
		a1 <- AR[1]
		a2 <- AR[2]
		b1 <- MA[1]
		VAR <- e_var*(a1*b1+a1*a2*b1+(1-a2)*(b1*(a1+b1)+1))/(1-a2-a1^2-a1^2*a2-a2^2*(1-a2))
		LAG1 <- (a1 + (1/VAR)*b1*e_var)/(1-a2)
		
	#ARMA(1,2)
	} else if (length(AR)==1 & length(MA)==2 & NonStationary==FALSE) {
		a1 <- AR[1]
		b1 <- MA[1]
		b2 <- MA[2]
		VAR <- e_var*(a1*b1+b2*a1^2+b2*b1*a1+1+b1*a1+b1^2+b2^2+b2*a1^2+b2*a1*b1)/(1-a1^2)
		LAG1 <- a1 + e_var*(b1+b2*a1+b2*b1)/VAR
		
	#ARMA(2,2)
	} else if (length(AR)==2 & length(MA)==2 & NonStationary==FALSE) {
		a1 <- AR[1]
		a2 <- AR[2]
		b1 <- MA[1]
		b2 <- MA[2]
		X <- 1-a2
		VAR <- e_var*(a1*b1 + a1^2*b2 + a1*b1*b2 + a1*a2*b1 + a1^2*a2*b2 + a1*a2*b1*b2 + X*(a2*b2 + 1 + b1*a1 + b1^2 + b2^2 + a2*b2 + b2*a1^2 + b2*a1*b1)) / (1 - a2 - a1^2 - a2*a1^2 - a2^2*(1-a2))
		LAG1 <- (a1*VAR + b1*e_var + b2*a1*e_var + b2*b1*e_var)/(VAR*(1-a2))
		
	} else if (NonStationary==FALSE) {
		winDialog(type="ok",message="We only accept models with order less than or equal to ARMA(2,2)")
		
	} else {
		VAR <- 9999999999999999999 
		LAG1 <- 999999999999999999
	}
	
	
	return(list(VAR,LAG1,e_var,a1,a2,b1,b2,NonStationary))
}

VAR_COV_ARMA_Alteration <- function(parameters,num_AR,num_MA,VAR_TARGET,LAG1_TARGET,e_var_org,AR_org,MA_org,Other_e_var,Other_AR,Other_MA,cur_dir) {
	#get variance and lag 1 autocorrelation of noise model
	RESULT <- VAR_COV_ARMA(parameters,num_AR,num_MA)
	VAR <- RESULT[[1]]
	LAG1 <- RESULT[[2]]
	e_var <-RESULT[[3]] 
	a1 <- RESULT[[4]]
	a2 <- RESULT[[5]]
	b1 <- RESULT[[6]]
	b2 <- RESULT[[7]]
	NonStationary <- RESULT[[8]]
	FINAL_VAR <- VAR
	FINAL_LAG1 <- LAG1
	
	#get variance and lag 1 autocorrelation of additional models if they exist
	num_low_freq_models <- length(Other_e_var)
	if (num_low_freq_models>0) {
		VAR_OTHER <- array(NA,num_low_freq_models)
		LAG1_OTHER <- array(NA,num_low_freq_models)
		LAG1_OTHER_WEIGHTED <- array(NA,num_low_freq_models)
		for (j in 1:num_low_freq_models) {
			cur_AR <- Other_AR[[j]]
			cur_MA <- Other_MA[[j]]
			cur_parameters <- c(Other_e_var[[j]],cur_AR,cur_MA)
			cur_num_AR <- length(cur_AR)
			cur_num_MA <- length(cur_MA)
			cur_result <- VAR_COV_ARMA(cur_parameters,cur_num_AR,cur_num_MA)
			VAR_OTHER[j] <- cur_result[[1]]
			LAG1_OTHER[j] <- cur_result[[2]]
			LAG1_OTHER_WEIGHTED[j] <- cur_result[[2]]*cur_result[[1]]		
		}
		FINAL_VAR <- sum(VAR_OTHER,VAR)
		FINAL_LAG1 <- sum(LAG1_OTHER_WEIGHTED,(LAG1*VAR))/FINAL_VAR		
	}
	
	e_dif <- e_var - e_var_org
	A1_dif <- 0 ; A2_dif <- 0; B1_dif <- 0; B2_dif <- 0
	if (num_AR==1 & NonStationary==FALSE) {
    A1_dif <- a1 - AR_org[1]
  }
	if (num_AR==2 & NonStationary==FALSE) {
    A1_dif <- a1 - AR_org[1]
    A2_dif <- a2 - AR_org[2]
  }
	if (num_MA==1 & NonStationary==FALSE) {
    B1_dif <- b1 - MA_org[1]
  }
	if (num_MA==2 & NonStationary==FALSE) {
    B1_dif <- b1 - MA_org[1]
    B2_dif <- b2 - MA_org[2]
  }
	Penalty <- abs(e_dif) + abs(A1_dif) + 10*abs(B1_dif) + 100*abs(A2_dif) + 1000*abs(B2_dif)
	
	distance <- 100000000*sqrt(((VAR_TARGET - FINAL_VAR)/VAR_TARGET)^2 + 
                               ((LAG1_TARGET-FINAL_LAG1)/LAG1_TARGET)^2)

	return(distance + Penalty)
}

Adjust_WARM_Model <- function (Stdev_Change,Lag1_Change,e_var_org,AR_org,MA_org,Other_e_var,Other_AR,Other_MA,cur_dir) {
  #e_var_org: the variance for the noise in the ARMA model to be adjusted
	#AR_org: a vector of AR coefficients for the ARMA model to be adjusted
	#MA_org: a vector of MA coefficients for the ARMA model to be adjusted
	#Other_e_var: a list of the variances for the noise in low frequency ARMA models. These are not adjusted. 
	#Other_AR: AR components of low frequency ARMA models. This is a list of vectors. These are not adjusted.
	#Other_MA: MA components of low frequency ARMA models. This is a list of vectors. These are not adjusted.
	#cur_dir: directory with adjustment functions
	
# 	setwd(cur_dir)
# 	source("Variance_Autocorrelation_ARMA_Alteration.r")
# 	source("Variance_Autocorrelation_ARMA.r")
# 	source("SCE_ALGORITHM.r")

	num_AR <- length(AR_org)
	num_MA <- length(MA_org)
	e_var_min <- 0.001*e_var_org
	e_var_max <- 100*e_var_org
	parameters <- c(e_var_org,AR_org,MA_org)
	RESULT <- VAR_COV_ARMA(parameters,num_AR,num_MA)
	VAR <- RESULT[[1]]
	LAG1 <- RESULT[[2]]
	FINAL_VAR <- VAR
	FINAL_LAG1 <- LAG1
	
	#get variance and lag 1 autocorrelation of additional models if they exist
	num_low_freq_models <- length(Other_e_var)
	VAR_OTHER <- 0
	LAG1_OTHER_WEIGHTED <- 0
	if (num_low_freq_models>0) {
		VAR_OTHER <- array(NA,num_low_freq_models)
		LAG1_OTHER <- array(NA,num_low_freq_models)
		LAG1_OTHER_WEIGHTED <- array(NA,num_low_freq_models)
		for (j in 1:num_low_freq_models) {
			cur_AR <- Other_AR[[j]]
			cur_MA <- Other_MA[[j]]
			cur_parameters <- c(Other_e_var[[j]],cur_AR,cur_MA)
			cur_num_AR <- length(cur_AR)
			cur_num_MA <- length(cur_MA)
			cur_result <- VAR_COV_ARMA(cur_parameters,cur_num_AR,cur_num_MA)
			VAR_OTHER[j] <- cur_result[[1]]
			LAG1_OTHER[j] <- cur_result[[2]]
			LAG1_OTHER_WEIGHTED[j] <- cur_result[[2]]*cur_result[[1]]		
		}
		FINAL_VAR <- sum(VAR_OTHER,VAR)
		FINAL_LAG1 <- sum(LAG1_OTHER_WEIGHTED,(LAG1*VAR))/FINAL_VAR		
	}

	VAR_TARGET <- (Stdev_Change*sqrt(FINAL_VAR))^2
	LAG1_TARGET <- max(min(Lag1_Change*FINAL_LAG1,.99),-.99)

	if (length(AR_org)==0 & length(MA_org)==0) {
		start_par <- c(e_var_org)	#e_var, a1
		lowerb <- c(e_var_min)
		upperb <- c(e_var_max)
	} else if (length(AR_org)==1 & length(MA_org)==0) {
		start_par <- c(e_var_org, AR_org[1])	#e_var, a1
		lowerb <- c(e_var_min,-1)
		upperb <- c(e_var_max,1)
	} else if (length(AR_org)==2 & length(MA_org)==0) {
		start_par <- c(e_var_org, AR_org[1], AR_org[2])	 #e_var, a1, a2
		lowerb <- c(e_var_min,-1,-1)
		upperb <- c(e_var_max,1,1)	
	} else if (length(AR_org)==0 & length(MA_org)==1) {
		start_par <- c(e_var_org,MA_org[1])	#e_var, b1
		lowerb <- c(e_var_min,-1)
		upperb <- c(e_var_max,1)  	
	} else if (length(AR_org)==0 & length(MA_org)==2) {
		start_par <- c(e_var_org,MA_org[1],MA_org[2])	#e_var, b1, b2
		lowerb <- c(e_var_min,-1,-1)
		upperb <- c(e_var_max,1,1)  	
	} else if (length(AR_org)==1 & length(MA_org)==1) {
		start_par <- c(e_var_org,AR_org[1],MA_org[1])	#e_var, a1, b1
		lowerb <- c(e_var_min,-1,-1)
		upperb <- c(e_var_max,1,1)  	
	} else if (length(AR_org)==2 & length(MA_org)==1) {
		start_par <- c(e_var_org,AR_org[1],AR_org[2],MA_org[1])	#e_var, a1, a2, b1
		lowerb <- c(e_var_min,-1,-1,-1)
		upperb <- c(e_var_max,1,1,1)	
	} else if (length(AR_org)==1 & length(MA_org)==2) {	
		start_par <- c(e_var_org,AR_org[1],MA_org[1],MA_org[2])	#e_var, a1, b1, b2
		lowerb <- c(e_var_min,-1,-1,-1)
		upperb <- c(e_var_max,1,1,1)
	} else if (length(AR_org)==2 & length(MA_org)==2) {
		start_par <- c(e_var_org,AR_org[1],AR_org[2],MA_org[1],MA_org[2])	#e_var, a1, a2, b1, b2
		lowerb <- c(e_var_min,-1,-1,-1,-1)
		upperb <- c(e_var_max,1,1,1,1)
	} else {
		winDialog(type="ok",message="We only accept models with order less than or equal to ARMA(2,2)")	
	}


	result <- SCEoptim(VAR_COV_ARMA_Alteration,start_par,num_AR,num_MA,VAR_TARGET,LAG1_TARGET,e_var_org,AR_org,MA_org,Other_e_var,Other_AR, Other_MA,cur_dir,lower=lowerb,upper=upperb)
	new_parameters <- result$par
	RESULT <- VAR_COV_ARMA(new_parameters,num_AR,num_MA)
	VAR_UPDATE <- sum(VAR_OTHER,RESULT[[1]])
	LAG1_UPDATE <- sum(LAG1_OTHER_WEIGHTED,(RESULT[[2]]*RESULT[[1]]))/VAR_UPDATE		
		
	Percent_Deviation_Var <- round((VAR_TARGET - VAR_UPDATE)*100/VAR_TARGET,2)
	Percent_Deviation_Lag1 <- round((LAG1_TARGET - LAG1_UPDATE)*100/VAR_TARGET,2)
	if (Percent_Deviation_Var>5) {print("Warning: Adjusted annual variance is >5% different than target annual variance")}
	if (Percent_Deviation_Lag1>5) {print("Warning: Adjusted annual lag-1 autocorrelation is >5% different than target annual lag-1 autocorrelation")}


	e_var_final <- new_parameters[1]
	AR_final <-NULL
	MA_final <- NULL
	if (num_AR>0) {AR_final <- new_parameters[2:(2+num_AR-1)]}
	if (num_MA>0) {MA_final <- new_parameters[(num_AR+2):(1+num_AR+num_MA)]}
	
  if (num_AR>0 & num_MA==0) {FINAL <- list(e_var_final,AR_final)}
  if (num_AR==0 & num_MA>0) {FINAL <- list(e_var_final,MA_final)}
  if (num_AR>0 & num_MA>0) {FINAL <- list(e_var_final,AR_final,MA_final)}
  if (num_AR==0 & num_MA==0) {FINAL <- list(e_var_final)}  
  
  return(FINAL)
}

## ----adjust noise model function-----------------------------------------
adjust.noise.model <- function(models, lag1.factor, sd.factor) {
  # lag1.factor - multiplicative adjustment factor for lag1 arma coefficients
  # sd.factor - multiplicative adjustment factor for noise standard deviation
  # models - list of arma models generated by WAVELET_DECOMPOSITION()
  
  # extract noise model and coefficients
  noise <- models[[1]]
  noise.e_var<- noise$sigma2
  noise.ar <- coef(noise)[grepl('ar', names(coef(noise)))]
  noise.ma <- coef(noise)[grepl('ma', names(coef(noise)))]
  
  # extract component models and coefficients
  comps.e_var <- list()
  comps.ar <- list()
  comps.ma <- list()
  for (l in 2:dim(wave.comp)[2]) {
    comp <- models[[l]]
  	comps.e_var[[l-1]] <- comp$sigma2
    comps.ar[[l-1]] <- coef(comp)[grepl('ar', names(coef(comp)))]
    comps.ma[[l-1]] <- coef(comp)[grepl('ma', names(coef(comp)))]
  }
  
  new_parameters <- Adjust_WARM_Model(sd.factor, lag1.factor, 
                                      noise.e_var, noise.ar, noise.ma,
                                      comps.e_var, comps.ar, comps.ma)
  
  noise.model.adj <- models[['NOISE']]
  noise.model.adj$sigma2 <- new_parameters[[1]]
  non_intercept <- which(names(models[['NOISE']]$coef)!="intercept")
  if (length(non_intercept)>0) {
  	for (cur_par in 1:length(non_intercept)) {
      noise.model.adj$coef[[cur_par]] <- new_parameters[[cur_par+1]]
    }
  }
  models[['NOISE']] <- noise.model.adj
  return(models)
}
adjust.noise.model(models=models, lag1.factor=1, sd.factor=1)

## ----show origin models--------------------------------------------------
models

## ----plot compare arma adjustment factors, fig.height=8------------------
par(mfrow=c(2,2))
adjust.noise.model(models=models, lag1.factor=1, sd.factor=0.5)[['NOISE']] %>%
  forecast(h=40) %>%
  plot(main="Lag1=1, sd=0.5")
adjust.noise.model(models=models, lag1.factor=1, sd.factor=1)[['NOISE']] %>%
  forecast(h=40) %>%
  plot(main="Lag1=1, sd=1")
adjust.noise.model(models=models, lag1.factor=0.5, sd.factor=0.5)[['NOISE']] %>%
  forecast(h=40) %>%
  plot(main="Lag1=2, sd=0.5")
adjust.noise.model(models=models, lag1.factor=0.5, sd.factor=1)[['NOISE']] %>%
  forecast(h=40) %>%
  plot(main="Lag1=2, sd=1")

## ----session-------------------------------------------------------------
sessionInfo()

