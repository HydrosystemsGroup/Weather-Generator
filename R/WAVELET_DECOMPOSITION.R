#' Wavelet decomposition
#'
#' @param x annual timeseries of observed climate variable (e.g. precipitation)
#' @param n.periods number of final periods as integer
#' @param sig.periods integer array of significance periods
#' @param n.comp.periods number of comparison periods
#' @param plot_flag boolean flag to show plots
#' @export
#' @examples
#' wavelet_decomposition(x=x, n.periods=1, sig.periods=seq(8,15), n.comp.periods=8, plot_flag=TRUE)

wavelet_decomposition <- function(x, n.periods, sig.periods, n.comp.periods, plot_flag=TRUE) {

  # wavelet analysis ----

  #....construct time series to analyze, pad if necessary
  CURRENT_CLIMATE_VARIABLE_org <- x
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

  # significance testing ----
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

  LOW_FREQUENCY_COMPONENTS <- array(0,c(length(x),n.periods))
  for (i in 1:n.periods) {
    CUR_PERIODS <- sig.periods[1:n.comp.periods[i]]
    if (i>1) {CUR_PERIODS <- sig.periods[(1 + (i-1)*n.comp.periods[i-1]):(n.comp.periods[i] + (i-1)*n.comp.periods[i-1])]}
    sj <- scale[CUR_PERIODS]
    #for Morlet Wavelet with freq = 6
    Cdelta <- .776
    w0_0 <- pi^(-1/4)
    if (length(CUR_PERIODS)>1) {LOW_FREQUENCY_COMPONENTS[,i] <- apply(sd(CURRENT_CLIMATE_VARIABLE_org)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[CUR_PERIODS,]/sqrt(sj),FUN=sum,c(2))}
    if (length(CUR_PERIODS)==1) {LOW_FREQUENCY_COMPONENTS[,i] <- sd(CURRENT_CLIMATE_VARIABLE_org)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[CUR_PERIODS,]/sqrt(sj)}
  }

  NOISE <- CURRENT_CLIMATE_VARIABLE_org - apply(LOW_FREQUENCY_COMPONENTS,FUN=sum,c(1))

  if (plot_flag) {
    par(mfrow=c((2+n.periods),1),font.axis=2,font.lab=2)
    plot(x,type="l",main="ORIGINAL DATA",xlab="TIME (YEARS)",ylab="")
    for (k in 1:n.periods){
      plot(LOW_FREQUENCY_COMPONENTS[,1],type="l",main=paste("COMPONENT",k),xlab="TIME (YEARS)",ylab="")
    }
    plot(NOISE,type="l",main="NOISE",xlab="TIME (YEARS)",ylab="")
  }
  return(cbind(NOISE,LOW_FREQUENCY_COMPONENTS))
}

