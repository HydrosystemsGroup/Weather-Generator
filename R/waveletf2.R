#' Wavelet function that returns fourier_factor, cone of influence and degrees of freedom
#'
#' @param k wavelet number array
#' @param s scale parameter
#' @export
#' @examples
#' waveletf2(k,s)

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
  dofmin <- 2    # Degrees of freedom
  return(c(fourier_factor,coi,dofmin))
}
