#' @title EM Algorithm for Gaussian Mixture Model Combined with Measurement Error Model and Linear Model
#' @param Zm Surrogate data in main study
#' @param Xv True observations in validation study
#' @param Zv Surrogate data in validation study
#' @param y Responses in outcome model
#' @param lambda A vector of size K containing initial value of mixting probabilities
#' @param muK A matrix of size K * p containing initial values of component mean, K-means center is specified if NULL
#' @param sigK An array with size p*p*K containing K p*p matrix, each of which is initial values of component variance,
#' K-means variance is specified if NULL
#' @param alpha A vector of size p consisting of intercept in measurement error model
#' @param A A matrix of size p*p consisting of slope in measurement error model
#' @param sigE A matrix of size p*P consisting of variance of noise in measurement error model
#' @param beta A vector size K containing initial linear coefficients
#' @param sigma2 A scalar of initial variance of random noise in linear model
#' @param profile A vector containing profile parameters
#' @param tol The convergence criterion
#' @param maxit The maximum number of iterations
#' @param verb If true, then various updates are printed during each iteration of the algorithm
#' @param is.profile Used for computing profile variance. If true, beta and sigma2 will not be estimated
#' @return A list with the following elements:
#' \itemize{
#' \item{pi} The final mixing probabilities
#' \item{mu} The final mean vectors
#' \item{sigK} The final variance matrix of each component
#' \item{alpha} The final \code{alpha} in measurement error model
#' \item{A} The final \code{A} in measurement error model
#' \item{sigK} The final \code{sigK} in measurement error model
#' \item{beta} The final linear coefficients
#' \item{sigma2} The final variance of random noise
#' \item{posterior} Membership probability
#' \item{iter} Number of iterations
#' }
#' @export
#' 
EM.mixme = function (Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2, tol=1e-5, maxit=5000, verb=FALSE, is.profile=FALSE) {
  K <- length(lambda); nm <- nrow(Zm); nv <- nrow(Zv); p <- ncol(Zm)
  ## E(Zi|Gi=k), return a K*p matrix
  Ezg = function(muK, alpha, A) t(A %*% t(muK) + alpha)
  ## Var(Zi|Gi=k), return a p*p*K array
  Vzg = function(sigK, alpha, sigE) array(apply(sigK, 3, function(x) (A %*% x %*% t(A) + sigE)), dim=c(p,p,K))
  ## Pr(Gi=k|Zi,Yi), return a n*K matrix 
  tildeomg = function(Z, lambda, Ezg1, Vzg1, tildemuK1, tildesigK1, y, beta, sigma2) {
    # each cluster
    temp <- do.call(cbind, lapply(1:K, function(k) lambda[k] * mvtnorm::dmvnorm(Z, Ezg1[k,], Vzg1[,,k]) * dnorm(y, beta[k], sqrt(sigma2))))
    return(temp/rowSums(temp))
  }
  ## tildesigK, a p*p*K array
  ## tildemuK, a n*p*K array
  tildemuKsigK = function(Z, muK, sigK, alpha, A, sigE) {
    temp <- lapply(1:K, function(k) sigK[,,k] %*% t(A) %*% solve(A %*% sigK[,,k] %*% t(A) + sigE))
    tildesigK1 <- list2array(lapply(1:K, function(k) sigK[,,k] - temp[[k]] %*% A %*% sigK[,,k]))
    tildemuK1 <- list2array(lapply(1:K, function(k) t(muK[k,] + temp[[k]] %*% (t(Z) - alpha - c(A %*% as.matrix(muK[k,]))) )))
    return(list(tildemuK=tildemuK1, tildesigK=tildesigK1))
  }
  ## E(I(Gi=k)Xi|Zi,Yi), return a n*p matrix 
  Egxz = function(k, tildemuK1, tildeomg1, muk=0) t(t(tildemuK1[,,k]) - muk) * tildeomg1[,k]
  ## E(I(Gi=k)XiXi^T|Zi,Yi), return a p*p*n array 
  Egxxz = function(k, tildemuK1, tildesigK1, tildeomg1, muk=0)
    list2array(lapply(1:dim(tildemuK1)[1], function(i) tildeomg1[i,k]*(tildesigK1[,,k] + outer(tildemuK1[i,,k] - muk, tildemuK1[i,,k] - muk))))
  ## E(Xi|Zi,Yi), return a n*p matrix 
  Exz = function(tildemuK1, tildeomg1) Reduce("+", lapply(1:K, function(k) Egxz(k, tildemuK1, tildeomg1)))
  ## E(XiXi^T|Zi,Yi), return a p*p*n array 
  Exxz = function(tildemuK1, tildesigK1, tildeomg1) Reduce("+", lapply(1:K, function(k) Egxxz(k, tildemuK1, tildesigK1, tildeomg1)))
  #### M-step
  ## update A and alpha
  # hatalphaA = function(Xf, Zf, Exz1, Exxz1) {
  #   # Zi = tildeXi * Gamma + ei
  #   # Gamma = (alpha01,...alpha0p,A11...A1p,...Ap1,...App)
  #   ## generate a p*p^2 row-block diagonal matrix
  #   diag.block = function(x) kronecker(diag(rep(1,length(x))),matrix(x,nr=1))
  #   ## transform Xi to tildeXi
  #   Ex2tld = function(x) cbind(diag(rep(1,p)),diag.block(x))
  #   Exx2tld = function(x,xx) rbind(Ex2tld(x),
  #                               cbind(t(diag.block(x)),
  #                               matrix(Matrix::bdiag(replicate(p,xx,simplify = F)),nr=p^2)))
  #   Exxz1f = abind::abind(Exxz1,list2array(lapply(1:nrow(Xv), function(i) outer(Xv[i,],Xv[i,]))))
  #   t1 = lapply(1:nrow(Xf), function(i) Exx2tld(Xf[i,], Exxz1f[,,i]))
  #   t2 = lapply(1:nrow(Xf), function(i) t(Ex2tld(Xf[i,])) %*% c(Zf[i,]) )
  #   out = solve(Reduce("+",t1)) %*% Reduce("+",t2)
  #   return(list(alpha = out[1:p], A = matrix(out[-(1:p)],nr=p,byrow=TRUE)))
  # }
  # Another solution: linear equation system -- from matrix to vector
  # A_{p*p} h_p + alpha_p c_0 = b_p
  # A_{p*p} M_{p*p} + alpha_p k_p = D_{p*p}
  solve.mateq = function(h, c_0, b, M, k, D) {
    # h : p-vector # c_0 : constant # b: p-vector
    # M: p*p matrix # k: p-vector # D: p*p matrix
    # solve equation G x = e
    # G = (A11,A12,...A1p,...,Ap1,..App,alpha01,...alpha0p)
    p <- length(h)
    # generate a p*p^2 row-block diagonal matrix
    diag.block = function(x) kronecker(diag(rep(1, length(x))), matrix(x, nr=1))
    G <- cbind(diag.block(h), diag(c_0, p, p))
    for (j in 1:p) {
      temp <- cbind(diag.block(M[,j]), diag(k[j], p, p))
      G <- rbind(G, temp)
    }
    e <- c(b, c(D))
    out <- as.vector(solve(G, e))
    return(list(A=matrix(out[1:p^2], nrow=p, byrow=T), alpha=tail(out, p)))
  }
  hatalphaA = function(Xf, Zf, Exz1, Exxz1) {
    h <- colSums(Xf); c_0 <- nm + nv; b <- colSums(Zf)
    M <- Reduce("+", array2list(Exxz1)) + Reduce("+", lapply(1:nv, function(i) outer(Xv[i,], Xv[i,])))
    k <- colSums(Xf)
    D <- Reduce("+", lapply(1:(nm + nv), function(i) outer(Zf[i,], Xf[i,]) ))
    return(solve.mateq(h, c_0, b, M, k, D))
  }
  ## update sigE
  hatsigE = function(Xf, Zf, alpha, A, Exxz1) {
    t1 <- Reduce("+", array2list(outf(t(Zf) - alpha)))
    t2 <- Reduce("+", lapply(1:nrow(Zf), function(i) outer(as.vector(A %*% as.matrix(Xf[i,])), Zf[i,] - alpha) ))
    Exxz2 <- list2array( append(array2list(Exxz1), array2list(outf(t(Xv)))) )
    t3 <- Reduce("+", lapply(1:nrow(Zf), function(i) A %*% Exxz2[,,i] %*% t(A)))
    return((t1 - t2 - t(t2) + t3)/nrow(Zf))
  }
  ## update lambda
  hatlambda = function(tildeomg1) colSums(tildeomg1)/nm
  ## update muK
  hatmuK = function(tildemuK1, tildeomg1) do.call(rbind, lapply(1:K,function(k) colSums(Egxz(k,tildemuK1, tildeomg1))/sum(tildeomg1[,k])))
  ## update sigK
  hatsigK = function(tildemuK1, tildesigK1, tildeomg1, muK)
    list2array(lapply(1:K, function(k) Reduce("+",array2list(Egxxz(k,tildemuK1, tildesigK1, tildeomg1, muK[k,]))) / sum(tildeomg1[,k]) ))
  ## update beta
  hatbeta = function(tildeomg1, y) colSums(tildeomg1 * y) / colSums(tildeomg1)
  ##
  hatsigma2 = function(beta, tildeomg1) (sum(y^2) + sum(t(tildeomg1)*beta^2) - sum(2*y*colSums(t(tildeomg1)*beta))) /nm
  #### log-likelihood
  
  #### EM Iteration
  count <- 0
  diff <- 1
  loglik <- sum(obs_loglik.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2))
  while (count < maxit & diff > tol) {
    old_loglik <- loglik
    ## Update E Step
    Ezg1 <- Ezg(muK, alpha, A)
    Vzg1 <- Vzg(sigK, alpha, sigE)
    tildemuKsigK1 <- tildemuKsigK(Zm, muK, sigK, alpha, A, sigE)
    tildemuK1 <- tildemuKsigK1$tildemuK
    tildesigK1 <- tildemuKsigK1$tildesigK
    tildeomg1 <- tildeomg(Zm, lambda, Ezg1, Vzg1, tildemuK1, tildesigK1, y, beta, sigma2)
    Exz1 <- Exz(tildemuK1, tildeomg1)
    Exxz1 <- Exxz(tildemuK1, tildesigK1, tildeomg1)
    # combine Z and combine EX
    Zf <- rbind(Zm, Zv)
    Xf <- rbind(Exz1, Xv)
    ## Update M Step
    alphaA <- hatalphaA(Xf, Zf, Exz1, Exxz1)
    alpha <- alphaA$alpha
    A <- alphaA$A
    sigE <- hatsigE(Xf, Zf, alpha, A, Exxz1)
    lambda <- hatlambda(tildeomg1)
    muK <- hatmuK(tildemuK1, tildeomg1)
    sigK <- hatsigK(tildemuK1, tildesigK1, tildeomg1, muK)
    if (!is.profile) {
      # update estimate of profile parameters
      beta <- hatbeta(tildeomg1, y)
      sigma2 <- hatsigma2(beta, tildeomg1)
    }
    ## convergence criterion
    loglik <- sum(obs_loglik.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2))
    diff <- abs(loglik - old_loglik)
    count <- count + 1
    if(verb) cat("iteration = ", count, "Log-likelihood diff is ", diff, 
                 "Observed log-likelihood is ", loglik, "\n")
  }
  return(list(pi=lambda, mu=muK, sigK=sigK, alpha=alpha, A=A, sigE=sigE, beta=beta, sigma2=sigma2,
              posterior=tildeomg1, iter=count, loglik=loglik))
}

#' @title EM Algorithm for Gaussian Mixture Model Combined with Linear Model
#' @param X Surrogate data or observations in main study
#' @param y Responses in outcome model
#' @param lambda A vector of size K containing initial value of mixing probabilities
#' @param muK A matrix of size K * p containing initial values of component mean, K-means center is specified if NULL
#' @param sigK An array with size p*p*K containing K p*p matrix, each of which is initial values of component variance,
#' K-means variance is specified if NULL
#' @param beta A vector size K containing initial linear coefficients
#' @param sigma2 A scalar of initial variance of random noise in linear model
#' @param tol The convergence criterion
#' @param maxit The maximum number of iterations
#' @param verb If true, then various updates are printed during each iteration of the algorithm
#' @param is.profile Used for computing profile variance. If true, beta and sigma2 will not be estimated
#' @return A list with the following elements:
#' \itemize{
#' \item{pi} The final mixing probabilities
#' \item{mu} The final mean vectors
#' \item{sigK} The final variance matrix of each component
#' \item{beta} The final linear coefficients
#' \item{sigma2} The final variance of random noise
#' \item{posterior} Membership probability
#' \item{iter} Number of iterations
#' }
#' @export
#' 
EM.mix = function (X, y, lambda, muK, sigK, beta, sigma2, tol=1e-5, maxit=5000, verb=FALSE, is.profile=FALSE) {
  #### E-step
  K <- length(lambda); nm <- nrow(X)
  ## Pr(Gi=k|Yi) return a n*K matrix
  tildeomg = function(lambda, X, y, beta, sigma2, muK, sigK) {
    temp <- do.call(cbind, lapply(1:K, function(k)
      lambda[k] * mvtnorm::dmvnorm(X, muK[k,], sigK[,,k]) * dnorm(y, beta[k], sqrt(sigma2))))
    return(temp/rowSums(temp))
  }
  #### M-step
  ## update lambda
  hatlambda = function(tildeomg1) colSums(tildeomg1)/nm
  ## update muK
  hatmuK = function(tildeomg1) do.call(rbind, lapply(1:K, function(k) colSums(tildeomg1[,k] * X)/sum(tildeomg1[,k])))
  ## update sigK
  hatsigK = function(tildeomg1, muK)
    list2array(lapply(1:K, function(k) Reduce("+", array2list(outf(t(sqrt(tildeomg1[,k]) *t(t(X) - muK[k,]) )))) /sum(tildeomg1[,k]) ))
  ## update  beta
  hatbeta = function(tildeomg1, y) colSums(tildeomg1 * y) / colSums(tildeomg1)
  ##
  hatsigma2 = function(beta, tildeomg1) (sum(y^2) + sum(t(tildeomg1)*beta^2) - sum(2*y*colSums(t(tildeomg1)*beta))) /nm
  #### EM Iteration
  count <- 0
  diff <- 1
  loglik <- sum(obs_loglik.mix(X, y, lambda, muK, sigK, beta, sigma2))
  while (count < maxit & diff > tol){
    old_loglik <- loglik
    ## Update E Step
    tildeomg1 <- tildeomg(lambda, X, y, beta, sigma2, muK, sigK)
    ## Update M Step
    lambda <- hatlambda(tildeomg1)
    muK <- hatmuK(tildeomg1)
    sigK <- hatsigK(tildeomg1, muK)
    if (!is.profile) {
      # update estimate of profile parameters
      beta <- hatbeta(tildeomg1, y)
      sigma2 <- hatsigma2(beta, tildeomg1)
    }
    ## convergence criterion
    loglik <- sum(obs_loglik.mix(X, y, lambda, muK, sigK, beta, sigma2))
    diff <- abs(loglik - old_loglik)
    count <- count + 1
    if(verb) cat("iteration = ", count, "Log-likelihood diff is ", diff, 
                 "Observed log-likelihood is ", loglik, "\n")
  }
  return(list(pi=lambda, mu=muK, sigK=sigK, beta=beta, sigma2=sigma2, posterior=tildeomg1, iter=count, loglik=loglik))
}
