## ----load data, warning=FALSE, message=FALSE-----------------------------
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
ggplot(clim.wyr, aes(WYEAR, PRCP)) +
  geom_line() +
  labs(x="Water Year", y="Annual Precipitation (mm)")

## ----scott code, echo=FALSE----------------------------------------------
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

## ----scott---------------------------------------------------------------
wave.scott <- WAVELET_ANALYSIS(CLIMATE_VARIABLE=clim.wyr$PRCP, siglvl=0.90,
                               background_noise="white", plot_flag=TRUE)

## ----scott result--------------------------------------------------------
wave.scott

## ----annual wavelet with sig of 0.8--------------------------------------
WAVELET_ANALYSIS(CLIMATE_VARIABLE=clim.wyr$PRCP, siglvl=0.80,
                 background_noise="white", plot_flag=TRUE)

## ----bw------------------------------------------------------------------
library(biwavelet)

bw <- biwavelet::wt(d=cbind(1:length(clim.wyr$PRCP), clim.wyr$PRCP), 
                    dt=1, dj=1/4, max.scale=length(clim.wyr$PRCP), 
                    sig.level=0.90, lag1=0, sig.test=1)
str(bw)

## ----plot_bw-------------------------------------------------------------
par(mfrow=c(1,1))
plot(bw)

## ----bw gws--------------------------------------------------------------
bw$gws <- apply(bw$power, 1, mean)

## ----bw gws sig----------------------------------------------------------
bw$gws.sig <- wt.sig(d=clim.wyr$PRCP, dt=bw$dt, scale=bw$scale, sig.test=1, 
                     sig.level=0.90, dof=length(clim.wyr$PRCP)-bw$scale, 
                     mother='morlet', lag1=0.0)

## ----wavelet bw compare--------------------------------------------------
# compare periods
cbind(biwavelet=bw$period, scott=wave.scott[[3]])

# compare global spectrum (row-wise mean of power for biwavelet::wt())
cbind(biwavelet=bw$gws, scott=wave.scott[[1]])

# compare signif
cbind(biwavelet=bw$gws.sig$signif, scott=wave.scott[[2]])

## ----bw output-----------------------------------------------------------
par(mfrow=c(1,2))
plot(bw, plot.cb=TRUE, plot.phase=FALSE)
plot(bw$gws, bw$period, type="b",
     xlab="Global Wavelet Spectrum", ylab="Fourier Period (Years)",
     log="y", ylim=rev(range(bw$period)), xlim=range(c(bw$gws, bw$gws.sig$signif)))
lines(bw$gws.sig$signif, bw$period, lty=2, col='red')  	

## ----wavelet analysis plot-----------------------------------------------
WAVELET_ANALYSIS(CLIMATE_VARIABLE=clim.wyr$PRCP, siglvl=0.90,
                 background_noise="white", plot_flag=TRUE)

## ----wgen----------------------------------------------------------------
wave.wgen <- wavelet_analysis(x=clim.wyr$PRCP, sig.level=0.90, noise.type="white")
str(wave.wgen)

## ----plot wgen-----------------------------------------------------------
par(mfrow=c(1,2))
plot(wave.wgen, plot.cb=TRUE, plot.phase=FALSE)
plot(wave.wgen$gws, wave.wgen$period, type="b",
     xlab="Global Wavelet Spectrum", ylab="Fourier Period (Years)",
     log="y", ylim=rev(range(wave.wgen$period)), 
     xlim=range(c(wave.wgen$gws, wave.wgen$gws.sig$signif)))
lines(wave.wgen$gws.sig$signif, wave.wgen$period, lty=2, col='red')    

## ----session-------------------------------------------------------------
sessionInfo()

