#' Adjust Daily Values by Shifting Gamma Distribution
#'
#' @param x vector of values (e.g. daily precipitation)
#' @param months vector of months corresponding to x
#' @param mean_changes adjustment factor for changing mean (scalar or vector of length 12 for monthly)
#' @param cv_changes adjustment factor for changing coef. of variation (scalar or vector of length 12 for monthly)
#' @param exclude_zeros if TRUE (default), excludes zero values from adjustment
#' @export
adjust_daily_gamma <- function(x, months, mean_changes=1, cv_changes=1, exclude_zeros=TRUE) {
  # TODO: mean_change by state_prob_change

  stopifnot(length(x)==length(months))

  if (length(mean_changes) == 1) {
    mean_changes <- rep(mean_changes, 12)
  }
  stopifnot(length(mean_changes)==12)

  if (length(cv_changes) == 1) {
    cv_changes <- rep(cv_changes, 12)
  }
  stopifnot(length(cv_changes)==12)

  xnew <- x[]
  meta <- list()

  for (m in seq(1, 12)) {
    if (exclude_zeros) {
      idx <- which(x > 0 & months == m)
    } else {
      idx <- which(months == m)
    }

    x_orig <- x[idx]
    gamma_fit <- MASS::fitdistr(x_orig, 'gamma')$estimate
    gamma_adj <- c(shape=gamma_fit[['shape']]/(cv_changes[m])^2,
                   rate=gamma_fit[['rate']]/(mean_changes[m]*(cv_changes[m])^2))

    probs <- (1:length(idx))/(length(idx)+1)
    q_gamma_orig <- qgamma(probs, shape=gamma_fit[['shape']], rate=gamma_fit[['rate']])
    q_gamma_adj <- qgamma(probs, shape=gamma_adj[['shape']], rate=gamma_adj[['rate']])
    ratios <- q_gamma_adj/q_gamma_orig

    multipliers <- approx(x=q_gamma_orig, y=ratios, xout=x_orig, rule=2)$y
    x_adj <- x_orig*multipliers

    xnew[idx] <- x_adj

    meta[[m]] <- list(factors=c(mean=mean_changes[m], cv=cv_changes[m]),
                      distrib=rbind(fit=gamma_fit, adjust=gamma_adj),
                      exclude_zeros=exclude_zeros,
                      idx=idx,
                      orig=x_orig,
                      adj=x_adj,
                      ratios=multipliers)
  }

  list(monthly=meta,
       changes=rbind(mean=mean_changes, cv=cv_changes),
       original=x,
       adjusted=xnew)
}
