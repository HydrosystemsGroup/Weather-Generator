## ----load libraries, warning=FALSE, message=FALSE------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())
library(weatherGen)
library(forecast)

## ----load data-----------------------------------------------------------
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
ggplot(clim.wyr, aes(WYEAR, PRCP)) +
  geom_line() +
  labs(x="Water Year", y="Annual Precipitation (mm)")

## ----wavelet-------------------------------------------------------------
wt.wgen <- wavelet_analysis(x=clim.wyr$PRCP, sig.level=0.90, noise.type="white")

## ----plot signif, fig.width=4--------------------------------------------
plot(wt.wgen$gws, wt.wgen$period, type="b",
     xlab="Global Wavelet Spectrum", ylab="Fourier Period (Years)",
     log="y", ylim=rev(range(wt.wgen$period)), 
     xlim=range(c(wt.wgen$gws, wt.wgen$gws.sig$signif)))
lines(wt.wgen$gws.sig$signif, wt.wgen$period, lty=2, col='red') 

## ----signif periods------------------------------------------------------
signif.periods.idx <- which(wt.wgen$gws > wt.wgen$gws.sig$signif)
signif.periods <- wt.wgen$period[signif.periods.idx]
print(signif.periods)

## ----wgen comp-----------------------------------------------------------
comps.wgen <- wavelet_components(x=clim.wyr$PRCP, wt=wt.wgen, 
                                 n.periods=length(wt.wgen$period), 
                                 sig.periods=seq(1,length(wt.wgen$period)), 
                                 n.comp.periods=rep(1, length(wt.wgen$period)))
head(comps.wgen)

## ----plot_wgen_comps-----------------------------------------------------
df.comps.wgen <- as.data.frame(comps.wgen) %>%
  mutate(TIMESTEP=row_number()) %>%
  gather(COMPONENT, VALUE, -TIMESTEP) %>%
  filter(COMPONENT != "NOISE") %>%
  mutate(COMPONENT=str_replace_all(COMPONENT, 'COMPONENT ', ''),
         COMPONENT=as.numeric(as.character(COMPONENT)),
         IS.SIGNIF=COMPONENT %in% signif.periods.idx) %>%
  merge(data.frame(COMPONENT=seq(1, length(wt.wgen$period)), PERIOD=wt.wgen$period),
        by="COMPONENT", all.x=TRUE) %>%
  mutate(LABEL=paste0("Period = ", format(PERIOD, digits=3), ' yr'))

ggplot(df.comps.wgen, aes(TIMESTEP, VALUE, color=IS.SIGNIF)) +
  geom_line() +
  facet_wrap(~LABEL) +
  scale_color_manual('', 
                     values=c('TRUE'='red', 'FALSE'='black'), 
                     labels=c('TRUE'='Significant', 'FALSE'='Not Significant')) +
  theme(legend.position='bottom')

## ----plot_wgen_comps_combine---------------------------------------------
comps.combine.wgen <- filter(df.comps.wgen, IS.SIGNIF) %>%
  group_by(TIMESTEP) %>%
  summarise(COMPONENT=sum(VALUE))
ggplot(comps.combine.wgen, aes(TIMESTEP, COMPONENT)) +
  geom_line()

## ----wgen comps noise----------------------------------------------------
comps.combine.wgen <- mutate(comps.combine.wgen,
                             PRCP=clim.wyr[['PRCP']],
                             NOISE=PRCP-COMPONENT)
comps.combine.wgen %>%
  gather(GROUP, VALUE, -TIMESTEP) %>%
  mutate(GROUP=ordered(as.character(GROUP), levels=c("PRCP", "COMPONENT", "NOISE"))) %>%
  ggplot(aes(TIMESTEP, VALUE)) +
  geom_line() +
  facet_wrap(~GROUP, ncol=1, scales='free_y')

## ----wgen comps alt combine----------------------------------------------
comps.wgen.sig <- wavelet_components(x=clim.wyr$PRCP, wt=wt.wgen, 
                                     n.periods=1, 
                                     sig.periods=signif.periods.idx, 
                                     n.comp.periods=length(signif.periods.idx))

## ----plot_wgen_comps_alt_combine-----------------------------------------
as.data.frame(comps.wgen.sig) %>%
  mutate(PRCP=clim.wyr[['PRCP']],
         TIMESTEP=row_number()) %>%
  gather(GROUP, VALUE, -TIMESTEP) %>%
  mutate(GROUP=ordered(as.character(GROUP), levels=c("PRCP", "COMPONENT 1", "NOISE"))) %>%
  ggplot(aes(TIMESTEP, VALUE)) +
  geom_line() +
  facet_wrap(~GROUP, ncol=1, scales='free_y')

## ----scott code, echo=FALSE----------------------------------------------
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

## ----comps scott---------------------------------------------------------
comps.scott <- WAVELET_DECOMPOSITION(CLIMATE_VARIABLE=clim.wyr$PRCP, 
                                     NUM_FINAL_PERIODS=1, 
                                     ALL_SIG_PERIODS=c(signif.periods.idx), 
                                     NUM_PERIODS_ALL_COMPS=length(signif.periods.idx),
                                     plot_flag=TRUE)

## ----compare comps-------------------------------------------------------
cbind(scott=comps.scott[, 2], weatherGen=comps.combine.wgen$COMPONENT)

## ----session-------------------------------------------------------------
sessionInfo()

