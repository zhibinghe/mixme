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
#' @param is.diff If true, the distributions of data in main study and validation study are different
#' @return A list with the following elements
#' \itemize{
#' \item{Zm} Surrogate data in the main study
#' \item{Xm} True observations in the main study
#' \item{Zv} Surrogate data in the validation study
#' \item{Xv} True observations in the validation study
#' \item{y} The outcome of linear model
#' \item{cluster.m} The true cluster labels of units in the main study
#' }
#' @export
#' 
gendat = function(n1, n2, lambda, muK, sigK, alpha, A, sigE, beta, sigma2, is.diff=FALSE) {
  K <- length(lambda)
  p <- length(alpha)
  lambda <- lambda/sum(lambda)
  # Generate true observation
  genX = function(size, is.diff) {
    # random cluster allocation
    indclus <- sample(1:K, size, replace=T, prob=lambda)
    X <- matrix(NA, nrow=size, ncol=p)
    for (k in 1:K) X[which(indclus==k),] <- mvtnorm::rmvnorm(sum(indclus==k), muK[k,], sigK[,,k])
    if (is.diff) X <- mvtnorm::rmvnorm(size, muK[-K,], sigK[,,-K]) # only valid for two cluster case
    return(list(X=X, ind.cluster = indclus))
  }
  # convert cluster index to matrix form with dimension nm * K 
  indexMatrix = function(indvec) {
    mat <- matrix(0, nrow=length(indvec), ncol=K) 
    mat[cbind(seq_along(indvec), indvec)] <- 1
    return(mat)
  }
  Xv <- genX(n2, is.diff=is.diff)$X
  Xm.data <- genX(n1, is.diff=FALSE)
  Xm <- Xm.data$X
  ym <- colSums(t(indexMatrix(Xm.data$ind.cluster)) * beta) + rnorm(n1, 0, sd=sqrt(sigma2))
  Zv <- alpha + A %*% t(Xv) + t(mvtnorm::rmvnorm(n2, mean=rep(0,nrow(sigE)), sigE))
  Zm <- alpha + A %*% t(Xm) + t(mvtnorm::rmvnorm(n1, mean=rep(0,nrow(sigE)), sigE))
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
mix.lm = function (X, y, K, lambda, muK, sigK, beta, sigma2, tol=1e-5, maxit=5000, verb=TRUE, variance=FALSE) {
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
  temp <- NULL
  if (variance) temp <- mix.variance(X, y, model$pi, model$mu, model$sigK, model$beta, model$sigma2, maxit=maxit)
  model$variance <- temp
  return(model)
}

#' @title Regression Calibration Model
#' @param Zm Surogate data in the main study
#' @param Xv True observations in the validation study
#' @param Zv Surogate data in the validation study
#' @param y  A vector consisting of outcome of linear model
#' @param K number of clusters
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @inheritParams mix.lm
#' @return output of linear outcome model
#' @export
#' 
mereg.lm = function(Zm, Xv, Zv, y, K, tol=1e-5, maxit=5000, verb=TRUE) {
  Zm <- as.matrix(Zm)
  Xv <- as.matrix(Xv)
  Zv <- as.matirx(Zv)
  fit <- lm(Xv ~ Zv)
  predXm <- t(t(Zm %*% coef(fit)[-1,]) + coef(fit)[1,])
  model <- mixtools::mvnormalmixEM(predXm, k=K, epsilon=tol, maxit=maxit, verb=verb)
  reg <- lm(y ~ -1 + ., data=as.data.frame(cbind(y, model$posterior)))
  return(reg)
}

#' @title Gaussian Mixture Model Combined with Measurement Error and Linear model
#' @description simultaneously estimate all the parameters
#' @inheritParams mereg.lm
#' @inheritParams mereg.lm
#' @inheritParams mereg.lm
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
                    tol=1e-5, maxit=5000, verb=TRUE, variance=FALSE) {
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
  model <- EM.mixme(Zm, Xv, Zv, y, lambda, muK, sigK, alpha, A, sigE, beta, sigma2, tol, maxit, verb)
  #### profile variance 
  temp <- NULL
  if (variance) temp <- mixme.variance(Zm, Xv, Zv, y, model$pi, model$mu, model$sigK, model$alpha, model$A, model$sigE,
                         model$beta, model$sigma2, maxit=maxit)
  model$variance <- temp
  return(model)
}

#' @title Component selection in Gaussian mixture model
#' @inheritParams EM.mix
#' @param M initial number of components
#' @param alpha tuning parameter controlling degree of penalty
#' @inheritParams EM.mix
#' @inheritParams EM.mix
#' @inheritParams EM.mix
#' @inheritParams EM.mix
#' @inheritParams EM.mix
#' @inheritParams EM.mix
#' @return A list with the following elements:
#' \itemize{
#' \item{pi} The final mixing probabilities, only containing selected nonzero lambda
#' \item{mu} The final mean vectors
#' \item{sigK} The final variance matrix of each component
#' \item{posterior} Membership probability
#' \item{iter} Number of iterations
#' \item{loglik} log-likelihood of observed data
#' \item{bic_loglik} bic criteria 
#' }
#' @export
#' 
gmm.sel = function(X, M, alpha, lambda, muK, sigK, maxit=5000, tol=1e-4, verb=FALSE) {
  ##
  K <- length(lambda); nm <- nrow(X); p = ncol(X)
  df <- 1 + 1.5* + p^2/2
  #### E-step
  ## Pr(Gi=k|Yi) return a n*K matrix
  tildeomg = function(lambda, X, muK, sigK) {
    temp <- do.call(cbind, lapply(1:length(lambda), function(k)
      lambda[k] * mvtnorm::dmvnorm(X, muK[k,], sigK[,,k]) ))
    return(temp/rowSums(temp))
  }
  obs.loglik = function(X, lambda) {
    ## 
    temp <- lapply(1:length(lambda), function(k) lambda[k] * mvtnorm::dmvnorm(X, muK[k,], sigK[,,k]))
    temp <- do.call(cbind, temp)
    return(log(rowSums(temp)))
  }
  #### M-step
  ## update lambda
  # hatlambda = function(tildeomg1) colSums(tildeomg1)/nm # standard
  hatlambda = function(tildeomg1, Mhat) {
    temp <- (colSums(tildeomg1)/nm - alpha*df) / (1 - Mhat*alpha*df)
    temp <- temp[abs(temp) > tol]
    return(temp)
  }
  ## update muK
  hatmuK = function(tildeomg1) do.call(rbind, lapply(1:ncol(tildeomg1), function(k) colSums(tildeomg1[,k] * X)/sum(tildeomg1[,k])))
  ## update sigK
  hatsigK = function(tildeomg1, muK)
    list2array(lapply(1:ncol(tildeomg1), function(k) Reduce("+", array2list(outf(t(sqrt(tildeomg1[,k]) *(X - muK[k,]))))) / sum(tildeomg1[,k])))
  #### EM Iteration
  count <- 0
  diff <- 1
  loglik <- sum(obs.loglik(X, lambda))
  tildeomg1 <- tildeomg(lambda, X, muK, sigK)
  while (count < maxit & diff > tol){
    #
    old_loglik <- loglik
    Mhat <- length(lambda)
    ## Update E Step
    lambda <- hatlambda(tildeomg1, Mhat)
    print(lambda)
    if (length(lambda) == 0 | any(is.na(lambda))) break
    tildeomg1 <- tildeomg(lambda, X, muK, sigK)
    ## Update M Step
    muK <- hatmuK(tildeomg1)
    sigK <- hatsigK(tildeomg1, muK)
    ## convergence criterion
    loglik <- sum(obs.loglik(X, lambda))
    diff <- abs(loglik - old_loglik)
    bic_loglik <- loglik - 0.5*length(lambda)*df*log(nm)
    count <- count + 1
    if(verb) cat("iteration = ", count, "Log-likelihood diff is ", diff, 
                 "Observed log-likelihood is ", loglik, "\n")
  }
  if (length(lambda) !=0) lambda <- lambda/sum(lambda)
  return(list(pi=lambda/sum(lambda), mu=muK, sigK=sigK, posterior=tildeomg1, iter=count, loglik=loglik, bic_loglik=bic_loglik))
}


