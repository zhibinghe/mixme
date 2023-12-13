#' @title Variance-covariance Estimate by Profile Likelihood Method
#' @description The variance of \code{beta} and \code{sigma} is of interest, others are taken as nuisance parameters
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @export
#' @return A variance-covariance matrix of (\code{beta}, \code{sigma2})
#' 
mixme.variance = function(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2) {
  ##
  nm <- nrow(Zm); nv <- nrow(Xv)
  log.theta <- obs_loglik.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2)
  log.profile <- matrix(NA, nr=nm+nv, nc=K+1)
  theta.profile <- c(beta, sigma2)
  # compute profile log-likelihood
  for(j in 1:(K+1)) {
    hn <- rep(0, K+1)
    hn[j] <- 1/sqrt(nm+nv)
    theta.profile.new <- theta.profile + hn
    est.nuisance <- EM.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta=theta.profile.new[1:K], 
                             sigma2=theta.profile.new[K+1], is.profile=TRUE, maxit=maxit)
    log.profile[,j] <- obs_loglik.mixme(Zm, Xv, Zv, y, est.nuisance$pi, est.nuisance$mu, est.nuisance$sigK,
                                        est.nuisance$alpha, est.nuisance$A, est.nuisance$sigE,
                                        theta.profile.new[1:K], theta.profile.new[K+1]) 
  }
  diff.profile <- (log.profile - log.theta) * sqrt(nm+nv)
  covariance <- solve(t(diff.profile) %*% diff.profile)
  colnames(covariance) <- rownames(covariance) <- c(paste0("beta", 1:K), "sigma2")
  return(covariance)
}

#' @title Variance-covariance Estimate by Profile Likelihood Method
#' @description The variance of \code{beta} and \code{sigma} is of interest, others are taken as nuisance parameters
#' @inheritParams mix.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @inheritParams mixme.lm
#' @export
#' @return A variance-covariance matrix of (\code{beta}, \code{sigma2})
#'
mix.variance = function(X, y, lambda, muK, sigK, beta, sigma2) {
  # compute log-likelihood 
  log.theta <- obs_loglik.mix(X, y, lambda, muK, sigK, beta, sigma2)
  # compute profile log-likelihood
  nm <- nrow(X)
  log.profile <- matrix(NA, nr=nm, nc=K+1)
  theta.profile <- c(beta, sigma2)
  for(j in 1:(K+1)) {
    hn <- rep(0, K+1)
    hn[j] <- 1/sqrt(nm)
    theta.profile.new <- theta.profile + hn
    est.nuisance <- EM.mix(X, y, lambda, muK, sigK, beta=theta.profile.new[1:K], 
                           sigma2=theta.profile.new[K+1], is.profile=TRUE, maxit=maxit)
    log.profile[,j] <- obs_loglik.mix(X, y, est.nuisance$pi, est.nuisance$mu, est.nuisance$sigK,
                                      theta.profile.new[1:K], theta.profile.new[K+1])
  }
  diff.profile <- (log.profile - log.theta) * sqrt(nm)
  covariance <- solve(t(diff.profile) %*% diff.profile)
  colnames(covariance) <- rownames(covariance) <- c(paste0("beta", 1:K), "sigma2")
  return(covariance)
}  

