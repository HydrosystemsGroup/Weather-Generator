#' Apply Morlet wavelet function to Fourier frequencies
#'
#' @param k numeric array of Fourier frequencies at which to compute wavelet
#' @param s wavelet scale parameter
#' @export
#' @return
#' List containing
#'
#' \itemize{
#'   \item \code{daughter}: a numeric vector of the wavelet function
#'   \item \code{fourier_factor}: ratio of Fourier period to scale
#'   \item \code{coi}: cone of influence
#'   \item \code{dofmin}: degrees of freedom for each point in wavelet power
#' }
#'
#' @references
#' http://paos.colorado.edu/research/wavelets/software.html
#' @examples
#' waveletf(k,s)

wavelet_morlet <- function(k,s) {
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
  return(list(daughter=daughter,
              fourier_factor=fourier_factor,
              coi=coi,
              dofmin=dofmin))
}
