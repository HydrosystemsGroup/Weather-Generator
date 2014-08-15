#' Wavelet analysis of annual timeseries
#'
#' @param x annual values of climate timeseries (e.g. precipitation) as numeric array
#' @param sig significance level
#' @param noise.type type of background noise ("white" or "red")
#' @param plot.flag boolean flag to show plots
#' @export
#' @examples
#' wavelet_annual(x=c(354.2, 234.4, 343), sig=0.90, noise.type="white", plot.flag=TRUE)

wavelet_annual <- function(x, sig=0.90, noise.type="white", plot.flag=TRUE) {
  #....construct time series to analyze, pad if necessary
	x.var <- var(x)
	x.n <- length(x)
	z <- scale(x)
	z.var <- var(z)
	base2 <- floor(log(x.n)/log(2) + 0.4999)   # power of 2 nearest to N
	z <- c(z,rep(0,(2^(base2+1)-x.n)))
	n <- length(z)

	#Determine parameters for Wavelet analysis
	dt <- 1
	dj <- 0.25
	s0 <- 2*dt
	J <- floor((1/dj)*log((x.n*dt/s0),base=2))

	#....construct SCALE array & empty PERIOD & WAVE arrays
	scale <- s0*2^((0:J)*dj)
	period <- scale
	wave <- array(as.complex(0),c(J+1,n))  # define the wavelet array

	#....construct wavenumber array used in transform [Eqn(5)]
	k <- c(1:floor(n/2))
	k <- k*((2.*pi)/(n*dt))
	k <- c(0,k,-rev(k[1:floor((n-1)/2)]))

	f <- fft(z, inverse=FALSE)        #fourier transform of standardized precipitation

	# loop through all scales and compute transform
	for (a1 in 1:(J+1)) {
		wav <- wavelet_morlet(k,scale[a1])
    daughter <- wav[['daughter']]
		fourier_factor <- wav[['fourier_factor']]
		coi <- wav[['coi']]
		dofmin <- wav[['dofmin']]
		wave[a1,] <- fft(f*daughter, inverse=TRUE)/n  # wavelet transform[Eqn(4)]
	}

	period <- fourier_factor*scale
	coi <- coi*dt*c((0.00001),1:((x.n+1)/2-1),rev((1:(x.n/2-1))),(0.00001))  # COI [Sec.3g]
	wave <- wave[,1:x.n]  # get rid of padding before returning
	POWER <- abs(wave)^2
	GWS <- x.var*apply(POWER, FUN=mean, 1) #Global Wavelet Spectrum

  # significance testing ----
  # get the appropriate parameters [see Table(2)]
	k0 <- 6
	empir <- c(2,0.776,2.32,0.60)
	dofmin <- empir[1]     # Degrees of freedom with no smoothing
	Cdelta <- empir[2]     # reconstruction factor
	gamma_fac <- empir[3]  # time-decorrelation factor
	dj0 <- empir[4]       # scale-decorrelation factor
	if (noise.type=="white") {lag1 <- 0}    #for red noise background, lag1 autocorrelation = 0.72, for white noise background, lag1 autocorrelation = 0
	if (noise.type=="red") {lag1 <- .72}

	freq <- dt / period   # normalized frequency
	fft_theor <- (1-lag1^2) / (1-2*lag1*cos(freq*2*pi)+lag1^2)  # [Eqn(16)]
	fft_theor <- fft_theor  # include time-series variance
	dof <- dofmin

	#ENTIRE POWER SPECTRUM
	chisquare <- qchisq(sig, dof)/dof
	signif <- fft_theor*chisquare   # [Eqn(18)]
	sig95 <- ((signif))%o%(array(1,x.n))  # expand signif --> (J+1)x(N) array
	sig95 <- POWER / sig95         # where ratio > 1, power is significant

	#TIME_AVERAGED (GLOBAL WAVELET SPECTRUM)
	dof <- x.n - scale
	if (length(dof) == 1) {dof <- array(0,(J+1))+dof}
	dof[which(dof < 1)] <- 1
	dof <- dofmin*sqrt(1 + (dof*dt/gamma_fac / scale)^2 )
	tt <- which(dof < dofmin)
	dof[tt] <- dofmin
	chisquare_GWS <- array(NA,(J+1))
	signif_GWS <- array(NA,(J+1))
	for (a1 in 1:(J+1)) {
		chisquare_GWS[a1] <- qchisq(sig,dof[a1])/dof[a1]
		signif_GWS[a1] <- fft_theor[a1]*x.var*chisquare_GWS[a1]
	}

  # plotting ----
	period_lower_limit <- 0
	sig_periods <- which(GWS>signif_GWS & period > period_lower_limit)
	if (plot.flag) {
		par(mfrow=c(1,2),font.axis=2,font.lab=2,cex.axis=1.1,cex.lab=1.1)
		xx <- 1:x.n
		yy <- period
		image(x=xx,y=yy,z=t(log(POWER,base=2)),xlab="Time (Years)",ylab="Fourier Period (Years)",ylim=rev(range(yy)),log="y",col=(heat.colors(12)))
		lines(xx,coi,lty=2)
		contour(x=xx,y=yy,z=t(sig95),add=TRUE,levels=c(-99,1),labels="")

		xmin <- (min(GWS,signif_GWS))
		xmax <- (max(GWS,signif_GWS))
		plot(GWS,yy,type="b",xlim=c(xmin,xmax),xlab="Global Wavelet Spectrum",ylab="Fourier Period (Years)",log="y",ylim=rev(range(yy)))
		lines(signif_GWS,yy,col="red",lty=2)
# 		if (noise.type=="white") {tt <- paste(sig*100,"% confidence \nlevel for white-noise\n spectrum",sep="")}
# 		if (noise.type=="red") {tt <- paste(sig*100,"% confidence \nlevel for red-noise\n spectrum",sep="")}
# 		legend(.2*max(GWS),.65*max(yy),tt,col="red",lty=2,box.lty=0,box.lwd=0,cex=.8)
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
		par(mfrow=c(1,1))
	}

	return(list(spectrum=GWS,signif=signif_GWS,period=period))
}



