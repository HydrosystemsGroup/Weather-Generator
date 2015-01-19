#' Compute State Equilibrium Vector of Transition Matrix
#'
#' @param m markov transition matrix
#' @export
#' @examples
#' transitions <- mc_fit(x=sample(c('d', 'w', 'e'), size=720, replace=TRUE, prob=c(0.5, 0.3, 0.2)), months=rep(rep(seq(1, 12), each=30), times=2))
#' mc_state_equilibrium(transitions[[1]])
#'
mc_state_equilibrium <- function(m) {
  eig <- eigen(t(m))$vectors[,1]
  eig <- eig/sum(eig)
  names(eig) <- rownames(m)
  eig
}
