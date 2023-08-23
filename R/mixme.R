# Transform a list to an array
list2array = function(X) array(unlist(X), dim=c(dim(X[[1]]), length(X)))
# Transform an array to a list
array2list = function(X) lapply(seq(dim(X)[3]), function(i) X[,,i])
# Outer product of each column of a matrix, return a nr*nr*nc array
outf = function(mat) array(apply(mat, 2, function(x) outer(x,x)), dim=c(nrow(mat), nrow(mat), ncol(mat)))

#' @title Simulate Data for Gaussian Mixture Model Combined with Measurement Error and Linear Model
#' @param n1 Number of units in the main study
#' @param n2 Number of units in the validation study
#' @param lambda A vector of size K consisting of mixing probabilities
#' @param muK A matrix of size K*p containing Gaussian component mean
#' @param sigK An array of size p*p*K, each element is the Gaussian component variance
#' @param alpha A vector of size p consisting of intercept in measurement error model
#' @param A A matrix of size p*p consisting of slope in measurement error model
#' @param sigE A matrix of size p*P consisting of variance of noise in measurement error model
#' @param beta A vector of size K consisting of linear coefficients in linear model
#' @param sigma2 A scalar of variance of random noise in linear model
#' @param Tn1 A scalar or a vector of size n1, repeat observations of units in Main study 
#' @param is.diff If true, the distributions of data in main study and validation study are different
#' @return A list with the following elements
#' \itemize{
#' \item{Zm} Surrogate data in the main study
#' \item{Xm} True observations in the main study
#' \item{Zv} Surrogate data in the validation study
#' \item{Xv} True observations in the validation study
#' \item{y} The outcome of linear model
#' \item{cluster.m} The True cluster labels of units in the main study
#' }
#' @export
#' 
gendat = function(n1, n2, lambda, muK, sigK, alpha, A, sigE, beta, sigma2, Tn1=1, is.diff=FALSE) {
  K <- length(lambda)
  p <- length(alpha)
  lambda <- lambda/sum(lambda)
  # Generate true observation
  genX = function(size, is.diff, w=1) {
    # random cluster allocation
    indclus <- sample(1:K, size, replace=T, prob=lambda)
    if (length(w) == 1) indclus <- rep(indclus, each=w)
    else indclus  = rep(indclus, times=w)
    X <- matrix(NA,length(indclus), p)
    for (k in 1:K) X[which(indclus==k),] <- mvtnorm::rmvnorm(sum(indclus==k), muK[k,], sigK[,,k])
    if (is.diff) X <- mvtnorm::rmvnorm(size, muK[-K,], sigK[,,-K])
    return(list(X=X, ind.cluster = indclus))
  }
  # convert cluster index to matrix form with dimension nm * K 
  indexMatrix = function(indvec) {
    mat <- matrix(0, nrow = length(indvec), ncol = K) 
    mat[cbind(seq_along(indvec), indvec)] <- 1
    return(mat)
  }
  Xv <- genX(n2,is.diff = is.diff)$X
  Xm.data <- genX(n1, is.diff=FALSE, w=Tn1)
  Xm <- Xm.data$X
  ym <- colSums(t(indexMatrix(Xm.data$ind.cluster)) * beta) + rnorm(nrow(Xm), 0, sd = sqrt(sigma2))
  Zv <- alpha + A %*% t(Xv) + t(mvtnorm::rmvnorm(nrow(Xv),mean = rep(0,nrow(sigE)), sigE))
  Zm <- alpha + A %*% t(Xm) + t(mvtnorm::rmvnorm(nrow(Xm),mean = rep(0,nrow(sigE)), sigE))
  return(list(Zm=t(Zm), Xm=Xm, Zv=t(Zv), Xv=Xv, y=ym, cluster.m=Xm.data$ind.cluster))
}

#' @title Compute Membership Probability for Given Observation 
#' @param x A vector or a matrix, if x is a matrix, each row is taken to be a quantile
#' @param lambda A vector of size K consisting of mixing probabilities
#' @param muK A vector of Gaussian mean
#' @param sigK A matrix of Gaussin variance matrix
#' @return A vector of size K consisting of membership probabilities
#' @export 
#' 
prob.member = function(x, lambda, muK, sigK) {
  K <- length(lambda) # number of clusters
  out <- do.call(cbind,lapply(1:K,function(k) lambda[k] * mvtnorm::dmvnorm(x, muK[k,], sigK[,,k])))
  return(out/rowSums(out))
}

#' @title Relabel the Clusters of Membership Probability
#' @description return matched Y corresponding to X
#' @param X A matrix of size n*K consisting of membership probabilities
#' @param Y A matrix of size n*K consisting of membership probabilities
#' @return Reordered Y by its columns, same cluster order as \code{X} 
#' 
match.c = function(X, Y) {
  if (ncol(X) != ncol(Y)) stop("dimensions of X and Y must be the same")
  ## permutation
  permutations = function(n) {
    if(n==1) return(matrix(1))
    else {
      sp <- permutations(n-1)
      p <- nrow(sp)
      A <- matrix(nrow=n*p, ncol=n)
      for(i in 1:n) A[(i-1)*p+1:p,] <- cbind(i, sp+(sp>=i))
      return(A)
    }
  }
  perm <- permutations(ncol(X))
  err <- apply(perm, 1, function(ind) mean(abs(X - Y[,ind])))
  out <- Y[, perm[which.min(err),]]
  return(out)
}

#' @title Initialization for Various EM Algorithms
#' @param X Data in main study
#' @param K Number of clusters
#' @param Xv Observations in validation study
#' @param Zv Surrogate data in validation study
#' @export
#' @return A list with the following elements:
#' \itemize{
#' \item{pi} Initial estimate of mixing probabilities 
#' \item{mu} A K*p matrix consisting of initial estimate of Gaussian mean 
#' \item{sigK} An K*p*p array consisting of initial estimate of Gaussian variance
#' \item{alpha} Initial estimate of intercept in measurement error model
#' \item{A} Initial estimate of slope in measurement error model
#' \item{sigE} Initial estimate of variance in measurement error model
#' }
#' 
mix.init = function (X=NULL, K=NULL, Xv=NULL, Zv=NULL) {
  ## GMM initial based on K-means
  if (is.null(X) | is.null(K)) lambda = muK = sigK = NULL
  else {
    kmeans_mod <- kmeans(X, K)
    lambda <- kmeans_mod$size/sum(kmeans_mod$size)
    muK <- kmeans_mod$centers
    sigK <- lapply(1:K, function(k) cov(X[which(kmeans_mod$cluster==k),]))
    sigK <- list2array(sigK)
  }
  ## ME initial
  if (is.null(Xv) | is.null(Zv)) alpha = A = sigE = NULL
  else {
    lm_mod <- lm(Zv ~ Xv)
    alpha <- lm_mod$coef[1,]
    A <- lm_mod$coef[-1,]
    sigE <- t(lm_mod$residuals) %*% lm_mod$residuals / (nrow(Xv) - ncol(Xv) -1)
  }
  return(list(pi=lambda, mu=muK, sigK=sigK, alpha=alpha, A=A, sigE=sigE))
}

## observed log-likelihood of mixme.lm model
obs_loglik.mixme = function (Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2) {
  ## part 1: main study
  temp <- lapply(1:K, function(k) lambda[k] * dnorm(y, beta[k], sqrt(sigma2)) * 
                   mvtnorm::dmvnorm(Zm, A %*% muK[k,] + alpha, A %*% sigK[,,k] %*% t(A) + sigE))
  temp <- do.call(cbind, temp)
  loglik1 <- log(rowSums(temp))
  ## part 2: validation study
  loglik2 <- log(sapply(1:nrow(Zv), function(i) mvtnorm::dmvnorm(Zv[i,], A %*% Xv[i,] + alpha, sigE)))
  return(c(loglik1, loglik2))
}
## observed log-likelihood of mix.lm model
obs_loglik.mix = function (X, y, lambda, muK, sigK, beta, sigma2) {
  ## 
  temp <- lapply(1:K, function(k) lambda[k] * dnorm(y, beta[k], sqrt(sigma2)) * 
                   mvtnorm::dmvnorm(X, muK[k,], sigK[,,k]))
  temp <- do.call(cbind, temp)
  return(log(rowSums(temp)))
}

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
#' @param is.profile If true, profile variance-covariance matrix will be computed
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
  K <- length(lambda); nm <- nrow(Zm); nv <- nrow(Zv)
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
  hatsigma2 = function(sigma0, beta, beta0, tildeomg1) sigma0^2 + sum(t(tildeomg1) * (beta - beta0)^2)/nm
  #### log-likelihood
  
  #### EM Iteration
  count <- 0
  diff <- 1
  loglik <- sum(obs_loglik.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2))
  while (count < maxit & diff > tol) {
    #
    beta0 <- beta; sigma0 <- sqrt(sigma2)
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
      sigma2 <- hatsigma2(sigma0, beta, beta0, tildeomg1)
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
#' @param is.profile If true, estimate nuisance parameters based on profile likelihood method
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
    list2array(lapply(1:K, function(k) Reduce("+", array2list(outf(t(sqrt(tildeomg1[,k]) *(X - muK[k,]))))) / sum(tildeomg1[,k])))
  ## update  beta
  hatbeta = function(tildeomg1, y) colSums(tildeomg1 * y) / colSums(tildeomg1)
  ##
  hatsigma2 = function(sigma0, beta, beta0, tildeomg1) sigma0^2 + sum(t(tildeomg1) * (beta - beta0)^2)/nm
  #### EM Iteration
  count <- 0
  diff <- 1
  loglik <- sum(obs_loglik.mix(X, y, lambda, muK, sigK, beta, sigma2))
  while (count < maxit & diff > tol){
    #
    old_loglik <- loglik
    beta0 <- beta; sigma0 <- sqrt(sigma2)
    ## Update E Step
    tildeomg1 <- tildeomg(lambda, X, y, beta, sigma2, muK, sigK)
    ## Update M Step
    lambda <- hatlambda(tildeomg1)
    muK <- hatmuK(tildeomg1)
    sigK <- hatsigK(tildeomg1, muK)
    if (!is.profile) {
      # update estimate of profile parameters
      beta <- hatbeta(tildeomg1, y)
      sigma2 <- hatsigma2(sigma0, beta, beta0, tildeomg1)
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

#' @title Gaussian Mixture Model Combined with Linear Model
#' @description Simultaneously estimate parameters in Gaussian mixture model combined with linear model
#' @param X A matrix consisting of the data
#' @param y A vector consisting of outcome of linear model
#' @param K Number of components
#' @param lambda A vector of size K containing initial value of mixting probabilities
#' @param muK A matrix of size K * p containing initial values of component mean, K-means center is specified if NULL
#' @param sigK An array with size p*p*K containing K p*p matrix, each of which is initial values of component variance,
#' K-means variance is specified if NULL
#' @param beta A vector size K containing initial linear coefficients
#' @param sigma2 A scalar of initial variance of random noise in linear model
#' @param tol The convergence criterion
#' @param maxit The maximum number of iterations
#' @param verb If true, then various updates are printed during each iteration of the algorithm
#' @param variance If true, variance matrix of beta is computed
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
mix.lm = function (X, y, K, lambda, muK, sigK, beta, sigma2, tol=1e-5, maxit=5000, verb=TRUE, variance=TRUE) {
  X <- as.matrix(X)
  p <- ncol(X)
  nm <- nrow(X)
  ## GMM initial based on K-means
  if (missing(lambda) | missing(muK) | missing(sigK)) {
    kmeans_mod <- mix.init(X, K)
    lambda <- kmeans_mod$pi
    muK <- kmeans_mod$mu
    sigK <- kmeans_mod$sigK
  }
  model <- EM.mix(X, y, lambda, muK, sigK, beta, sigma2, tol, maxit, verb)
  #### profile variance 
  if (variance) {
    temp <- mix.variance(X, y, model$pi, model$mu, model$sigK, model$beta, model$sigma2)
    model$variance <- temp
  }
  return(model)
}

#' @title Regression Calibration Model
#' @param Zm Surogate data in the main study
#' @param Xv True observations in the validation study
#' @param Zv Surogate data in the validation study
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @return output of mixtools::mvnormalmixEM
#' @export
#' 
mereg = function(Zm, Xv, Zv, lambda, muK, sigK, tol=1e-5, maxit=5000, verb=FALSE) {
  mod.dat <- as.data.frame(cbind(Xv, Zv))
  colnames(mod.dat) <- c(paste("X", 1:ncol(Xv), sep=""), paste("Z", 1:ncol(Zv), sep=""))
  Zm <- as.data.frame(Zm)
  colnames(Zm) <- paste("Z", 1:ncol(Zm), sep="")
  fit <- lm(cbind(X1, X2) ~ ., mod.dat)
  predXm <- suppressWarnings(predict(fit, newdata=Zm))
  mod <- mixtools::mvnormalmixEM(predXm, lambda=lambda, mu=lapply(1:nrow(muK), function(k) muK[k,]),
                                 sigma=lapply(1:nrow(muK), function(k) sigK[,,k]),
                                 k=nrow(muK), epsilon=tol, maxit=maxit, verb=verb)
  return(mod)
}

#' @title Gaussian Mixture Model Combined with Measurement Error and Linear model
#' @description simultaneously estimate all the parameters
#' @inheritParams mereg
#' @inheritParams mereg
#' @inheritParams mereg
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams EM.mixme
#' @inheritParams EM.mixme
#' @inheritParams EM.mixme
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @param variance If true, variance of \code{beta} and \code{sigma2} are estimated based on profile likelihood method
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
mixme.lm = function(Zm, Xv, Zv, y, K, lambda, muK, sigK, alpha, A, sigE, beta, sigma2,
                    tol=1e-5, maxit=5000, verb=TRUE, variance=TRUE) {
  #### mixme EM algorithm 
  Zm <- as.matrix(Zm); Zv <- as.matrix(Zv); Xv <- as.matrix(Xv)
  p <- ncol(Xv); nm <- nrow(Zm); nv <- nrow(Zv)
  ## GMM initial
  if (missing(lambda) | missing(muK) | missing(sigK)) {
    kmeans_mod <- mix.init(X=Xv, K=K)
    lambda <- kmeans_mod$pi
    muK <- kmeans_mod$mu
    sigK <- kmeans_mod$sigK
  }
  ## ME initial
  if (missing(alpha) | missing(A) | missing(sigE)) {
    kmeans_mod <- mix.init(Xv=Xv, Zv=Zv)
    alpha <- kmeans_mod$alpha
    A <- kmeans_mod$A
    sigE <- kmeans_mod$sigE
  }
  ## LM initial
  model <- EM.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2, 
                    tol=tol, maxit=maxit, verb=verb, is.profile=FALSE)
  #### profile variance 
  temp <- mixme.variance(Zm, Xv, Zv, y, model$pi, model$mu, model$sigK, model$alpha, model$A, model$sigE,
                         model$beta, model$sigma2)
  model$variance <- temp
  return(model)
}

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
                             sigma2=theta.profile.new[K+1], is.profile=TRUE)
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
                             sigma2=theta.profile.new[K+1], is.profile=TRUE)
    log.profile[,j] <- obs_loglik.mix(X, y, est.nuisance$pi, est.nuisance$mu, est.nuisance$sigK,
                                      theta.profile.new[1:K], theta.profile.new[K+1])
  }
  diff.profile <- (log.profile - log.theta) * sqrt(nm)
  covariance <- solve(t(diff.profile) %*% diff.profile)
  colnames(covariance) <- rownames(covariance) <- c(paste0("beta", 1:K), "sigma2")
  return(covariance)
}  

