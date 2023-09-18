#Moments estimates of the Dirichlet Distribution

#### In order not to call the gtools, just copy the simple functions here to be used
"dirichmomcopy" <- function(X){
    mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2)-((mean(X[,1]))^2))
    return(matrix(mom))
  }

"rdirichletcopy" <- function(n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x / as.vector(sm)
}

"ddirichletcopy" <- function(x, alpha) {
  dirichlet1 <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- (alpha - 1) * log(x)
    s <- ifelse(alpha == 1 & x == 0, -Inf, s)
    exp(sum(s) - logD)
  }
  # make sure x is a matrix
  if (!is.matrix(x)) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    } else {
      x <- t(x)
    }
  }
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), byrow = TRUE)
  }
  if (any(dim(x) != dim(alpha))) {
    stop("Mismatch between dimensions of 'x' and 'alpha'.")
  }
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) {
    pd[i] <- dirichlet1(x[i, ], alpha[i, ])
  }
  # Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
  pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
  pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
  pd
}


BerOMNPP_MCMC1 <- function(n0,y0,n,y,prior_gamma,prior_p,gamma_ind_prop,
                           gamma_ini,nsample,burnin,thin, adjust = FALSE){
  node = c(nsample/4, nsample/2)
  k = length(n0) # number of historical datasets

  LogPostBergamma = function(x){
    cumsumx = cumsum(x[1:k])
    gammay0 = sum(cumsumx*y0)
    gammanmy0 = sum(cumsumx*(n0-y0))
    out = lbeta(gammay0+y+prior_p[1],n-y+gammanmy0+prior_p[2])-
      lbeta(gammay0+prior_p[1],gammanmy0+prior_p[2])+sum((prior_gamma-1)*log(x))
    return(out)
  }
  if(is.null(gamma_ini)){
    gamma_ini = rep(1/(k+1), k+1)
  }
  gamma_cur = gamma_ini
  gamma_sample = matrix(0, nrow = nsample, ncol = k+1)
  counter = 0
  niter = nsample*thin + burnin

  for (i in 1:niter){
    ## Independent Dirichlet Proposal with adjustment at selected iterations
    ii = (i-burnin)/thin
    if(isTRUE(adjust)){
      if(ii %in% node){
        gamma_ind_prop = c(dirichmomcopy(gamma_sample[1:ii, ]))
      }
    }
    gamma_prop = rdirichletcopy(1, alpha = gamma_ind_prop)
    lpost_prop = min(LogPostBergamma(gamma_prop), 1e100)
    lpost_prop = max(lpost_prop, -1e100)                           ## Lazy way To Avoid Overflow
    lpost_cur = min(LogPostBergamma(gamma_cur), 1e100)
    lpost_cur = max(lpost_cur, -1e100)
    safelprior_cur = min(log(ddirichletcopy(gamma_cur, alpha = gamma_ind_prop)), 1e100)
    safelprior_cur = max(safelprior_cur, -1e100)
    safelprior_prop = min(log(ddirichletcopy(gamma_prop, alpha = gamma_ind_prop)), 1e100)
    safelprior_prop = max(safelprior_prop, -1e100)
    logr = min(0, (lpost_prop-lpost_cur+safelprior_cur-safelprior_prop))

    if(runif(1) <= exp(logr)){
      gamma_cur = gamma_prop
      counter = counter+1
    }
    if(i > burnin & (i-burnin) %% thin==0) {
      gamma_sample[ii, ] = gamma_cur
    }
  }
  gamma_minus = gamma_sample[, 1:k]
  delta_sample = t(apply(gamma_minus, 1, cumsum))
  # Generate p conditional on gamma; vectorized
  s1 = colSums(y0 * t(delta_sample)) + y + prior_p[1]
  s2 = colSums((n0-y0) * t(delta_sample)) + n - y + prior_p[2]
  p = rbeta(nsample, shape1 = s1, shape2 = s2)
  ans = list(acceptrate=counter/niter,p=p,delta=delta_sample)
  return(ans)
}
