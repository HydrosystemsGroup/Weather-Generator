#' Solve transition matrix
#'
#' @param p00 transition probability from 0 to 0
#' @param p01 transition probability from 0 to 1
#' @param p02 transition probability from 0 to 2
#' @param p10 transition probability from 1 to 0
#' @param p11 transition probability from 1 to 1
#' @param p12 transition probability from 1 to 2
#' @param p20 transition probability from 2 to 0
#' @param p21 transition probability from 2 to 1
#' @param p22 transition probability from 2 to 2
#' @export
#' @examples
#' GET_PI(p00,p01,p02,p10,p11,p12,p20,p21,p22)

GET_PI <- function(p00,p01,p02,p10,p11,p12,p20,p21,p22) {
  pi0 <- 1
  pi1 <- 1
  pi2 <- 1
  P <- matrix(c(p00,p01,p02,pi0,p10,p11,p12,pi1,p20,p21,p22,pi2),byrow=FALSE,nrow=4)
  P[1,1] <- P[1,1] - 1
  P[2,2] <- P[2,2] - 1
  P[3,3] <- P[3,3] - 1
  B <- c(0,0,0,1)
  PI <- solve(P[c(1,2,4),],B[c(1,2,4)])
  return(PI)
}
