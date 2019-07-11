#### Calculate Values of the normalizing constant on selected knots
#### Input: D_0, ln L(theta|D_0)
#### Input: The form of ln L(theta|D_0)^{delta}*pi(theta) - 
####        used to do sampling from the power prior with fixed delta

#### To calculate log C(delta) with given inputs:
#### Log Likelihood Function of D0, and theta 
#### Matrix of theta (M*n) dimension 
#### The corresponding delta (n dimension vector)

#### To calculate loglik(D0, theta) for the whole list
#### Input: D0 (Vector or Matrix)
#### Input: thetalist: list of the theta from L(theta|D_0)^{delta}*pi(theta)
#### Input: ntheta: number of parameters

#### Output: llik: Matrix of values for loglik of D0 given Theta, 
####               with M (Monte Carlo size) row and k (number of knots) column
#### Output: dknot: knots of delta 

loglikNormD0 <- function(D0, thetalist, ntheta = 2){
  if(ntheta!= length(thetalist)) stop("Number of free-parameters must match number of elements in the list")
  loglik <- dnorm(D0, mean = thetalist$thetamat1, sd = sqrt(thetalist$thetamat2), log = T)
  class(loglik) <- c( "npp")
  return(loglik)
}

loglikBerD0 <- function(D0, thetalist, ntheta = 1){
  if(ntheta!= length(thetalist)) stop("Number of free-parameters must match number of elements in the list")
  nsize <- length(D0)
  success <- sum(D0)
  loglik <- dbinom(x = success, size = nsize, prob = thetalist$thetamat, log = T)
  class(loglik) <- c( "npp")
  return(loglik)
}



logCknot <- function(deltaknot, llikf0){
  if(!inherits(llikf0, "npp")) stop("'llikf0' must be a function of the class npp")
  loglikmat <- llikf0
  loglikmean <- apply(loglikmat, 2, mean)
  dvec <- c(deltaknot[1], diff(deltaknot))
  lgC_knot <- cumsum(dvec*loglikmean)
  return(lgC_knot)
}


logCdelta <- function(delta, deltaknot, lCknot){
  if(delta %in% deltaknot){
  index <- match(delta, deltaknot)
  ans <- lCknot[index]
  }else{
  diff <- deltaknot-delta
  index <- which(diff > 0)[1]
  rate <- (lCknot[index] - lCknot[index-1])/(deltaknot[index]-deltaknot[index-1])
  ans <- lCknot[index-1]+rate*(delta-deltaknot[index-1])
  }
  return(ans)
}
