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
  x.var <- var(x)
  x.n <- length(x)
  z <- scale(x)
  z.var <- var(z)
  base2 <- floor(log(x.n)/log(2) + 0.4999)   # power of 2 nearest to N
  z.mirror <- c(z, rep(0, (2^(base2+1) - x.n)))
  n <- length(z.mirror)

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

  f <- fft(z.mirror, inverse=FALSE)        #fourier transform of standardized precipitation

  # loop through all scales and compute transform
  for (a1 in 1:(J+1)) {
    wavelet <- wavelet_morlet(k,scale[a1])
    daughter <- wavelet[['daughter']]
    fourier_factor <- wavelet[['fourier_factor']]
    coi <- wavelet[['coi']]
    dofmin <- wavelet[['dofmin']]
    wave[a1,] <- fft(f*daughter, inverse=TRUE)/n  # wavelet transform[Eqn(4)]
  }

  period <- fourier_factor*scale
  coi <- coi*dt*c((0.00001),1:((x.n+1)/2-1),rev((1:(x.n/2-1))),(0.00001))  # COI [Sec.3g]
  wave <- wave[,1:x.n]  # get rid of padding before returning
  wave.power <- abs(wave)^2
  wave.spec <- x.var*apply(wave.power,FUN=mean,c(1)) #Global Wavelet Spectrum

  # significance testing ----
  background_noise <- "white"   #can be red or white

  siglvl <- 0.95
  # get the appropriate parameters [see Table(2)]
  k0 <- 6
  empir <- c(2,0.776,2.32,0.60)
  dofmin <- empir[1]     # Degrees of freedom with no smoothing
  Cdelta <- empir[2]     # reconstruction factor
  gamma_fac <- empir[3]  # time-decorrelation factor
  dj0 <- empir[4]        # scale-decorrelation factor
  if (background_noise=="white") {lag1 <- 0} # for red noise background, lag1 autocorrelation = 0.72
  if (background_noise=="red") {lag1 <- .72} # for white noise background, lag1 autocorrelation = 0

  freq <- dt / period   # normalized frequency
  fft_theor <- (1-lag1^2) / (1-2*lag1*cos(freq*2*pi)+lag1^2)  # [Eqn(16)]
  fft_theor <- fft_theor  # include time-series variance
  dof <- dofmin

  #ENTIRE POWER SPECTRUM
  chisquare <- qchisq(siglvl,dof)/dof
  signif <- fft_theor*chisquare         # [Eqn(18)]
  sig95 <- ((signif))%o%(array(1,x.n))  # expand signif --> (J+1)x(N) array
  sig95 <- wave.power / sig95           # where ratio > 1, power is significant

  #TIME_AVERAGED (GLOBAL WAVELET SPECTRUM)
  dof <- x.n - scale
  if (length(dof) == 1) {dof <- array(0,(J+1))+dof}
  dof[which(dof < 1)] <- 1
  dof <- dofmin*sqrt(1 + (dof*dt/gamma_fac / scale)^2 )
  tt <- which(dof < dofmin)
  dof[tt] <- dofmin
  chisquare.spec <- array(NA,(J+1))
  signif.spec <- array(NA,(J+1))
  for (a1 in 1:(J+1)) {
    chisquare.spec[a1] <- qchisq(siglvl,dof[a1])/dof[a1]
    signif.spec[a1] <- fft_theor[a1]*x.var*chisquare.spec[a1]
  }

  # ------------------------------------
  freq.comps <- array(0,c(length(x),n.periods))
  for (i in 1:n.periods) {
    cur.periods <- sig.periods[1:n.comp.periods[i]]
    if (i>1) {
      cur.periods <- sig.periods[(1 + (i-1)*n.comp.periods[i-1]):(n.comp.periods[i] + (i-1)*n.comp.periods[i-1])]
    }
    sj <- scale[cur.periods]
    #for Morlet Wavelet with freq = 6
    Cdelta <- .776
    w0_0 <- pi^(-1/4)
    if (length(cur.periods)>1) {
      freq.comps[,i] <- apply(sd(x)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[cur.periods,]/sqrt(sj),FUN=sum,c(2))
    }
    if (length(cur.periods)==1) {
      freq.comps[,i] <- sd(x)*(dj*sqrt(dt)/(Cdelta*w0_0))*Re(wave)[cur.periods,]/sqrt(sj)
    }
  }

  noise <- x - apply(freq.comps,FUN=sum,c(1))

  if (plot_flag) {
    par(mfrow=c((2+n.periods),1))
    plot(x,type="l",main="ORIGINAL DATA",xlab="TIME (YEARS)",ylab="")
    for (k in 1:n.periods){
      plot(freq.comps[,1],type="l",main=paste("COMPONENT",k),xlab="TIME (YEARS)",ylab="")
    }
    plot(noise,type="l",main="NOISE",xlab="TIME (YEARS)",ylab="")
    par(mfrow=c(1,1))
  }

  colnames(freq.comps) <- paste("COMPONENT", 1:n.periods, sep=' ')

  return(cbind(NOISE=noise, freq.comps))
}

