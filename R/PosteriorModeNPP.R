
ModeDeltaBerNPP <- function(Data.Cur, Data.Hist,
                            CompStat = list(n0 = NULL, y0 = NULL, n1 = NULL, y1 = NULL), npoints = 1000,
                            prior = list(p.alpha = 1, p.beta = 1,
                                         delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    y0 <- Data.Hist[2]; n0 <- Data.Hist[1]
    y1 <- Data.Cur[2]; n1 <- Data.Cur[1]
  }else{
    y0 <- CompStat$y0
    n0 <- CompStat$n0
    y1 <- CompStat$y1
    n1 <- CompStat$n1
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen <- lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta)-
    lbeta(x*y0 + prior$p.alpha, x*(n0-y0) + prior$p.beta) +
    (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)

  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}


ModeDeltaNormalNPP <- function(Data.Cur, Data.Hist,
                             CompStat = list(n0 = NULL, mean0 = NULL, var0 = NULL, n1 = NULL, mean1 = NULL, var1 = NULL),
                             npoints = 1000,
                             prior = list(a = 1.5, delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    mean0 <- mean(Data.Hist)
    n0 <- length(Data.Hist)
    var0 <- var(Data.Hist)*(n0-1)/n0
    mean1 <- mean(Data.Cur)
    n1 <- length(Data.Cur)
    var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    mean0 <- CompStat$mean0
    n0 <- CompStat$n0
    var0 <- CompStat$var0
    mean1 <- CompStat$mean1
    n1 <- CompStat$n1
    var1 <- CompStat$var1
  }
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)

  propDen <-  (x*n0/2+prior$a+prior$delta.alpha-2)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma((x*n0+n1-3)/2+prior$a)-lgamma((x*n0-3)/2+prior$a)-
    ((x*n0+n1-3)/2+prior$a)*log((x*n1*(mean0-mean1)^2)/((x*n0+n1)*var0)+x+(n1*var1)/(n0*var0))

  deltamin = max(0, (3-2*prior$a)/n0)
  propDen = ifelse(x>deltamin, propDen, log(.Machine$double.xmin))

  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}




ModeDeltaPoisNPP <- function(Data.Cur, Data.Hist, CompStat = list(n0 = NULL, mean0 = NULL, n1 = NULL, mean1 = NULL),
                             npoints = 1000, prior = list(lambda.shape = 1/2, lambda.scale = 100,
                                                          delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    n0 <- length(Data.Hist)
    n1 <- length(Data.Cur)
    sum0 <- sum(Data.Hist)
    sum1 <- sum(Data.Cur)
  }else{
    n0 <- CompStat$n0
    n1 <- CompStat$n1
    sum0 <- n0*CompStat$mean0
    sum1 <- n1*CompStat$mean1
  }

  prior.lambda.rate <- 1/prior$lambda.scale
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)
  propDen <-  (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma(sum1+x*sum0+prior$lambda.shape)-
    lgamma(x*sum0+prior$lambda.shape)+
    (x*sum0+prior$lambda.shape)*log(n0*x+prior.lambda.rate)-
    (x*sum0+sum1+prior$lambda.shape)*log(n0*x+n1+prior.lambda.rate)

  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}




ModeDeltaMultinomialNPP <- function(Data.Cur, Data.Hist, CompStat = list(n0 = NULL, n1 = NULL),  npoints = 1000,
                                    prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5), delta.alpha = 1, delta.beta = 1))
{
  if(missing(CompStat)){
    if(length(Data.Cur) != length(Data.Cur)) stop("Number of Categories Must Be Equal!")
    n0 <- Data.Hist
    n1 <- Data.Cur
  }else{
    if(length(CompStat$n0) != length(CompStat$n1)) stop("Number of Categories Must Be Equal!")
    n0 <- CompStat$n0
    n1 <- CompStat$n1
  }

  n0sum <- sum(n0); n1sum <- sum(n1)
  alphasum <- sum(prior$theta.dir.alpha)
  x <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)

  s1 <- c(); s2 <- c()
  for(i in 1:length(x)){
    s1[i] <- sum(lgamma(n0*x[i]+n1+prior$theta.dir.alpha))
    s2[i] <- sum(lgamma(n0*x[i]+prior$theta.dir.alpha))
  }
  propDen <- ((prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
                lgamma(alphasum+x*n0sum)+s1-s2-
                lgamma(alphasum+x*n0sum+n1sum))

  NPPpostLik <- data.frame(x = x, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}




ModeDeltaLMNPP <- function(y.Cur, y.Hist, x.Cur = NULL, x.Hist = NULL, npoints = 1000,
                           prior = list(a = 1.5, b = 0, mu0 = 0, Rinv = matrix(1, nrow = 1),
                                        delta.alpha = 1, delta.beta = 1))
{
  a = prior$a; b = prior$b
  mu0 = as.matrix(prior$mu0); R = solve(prior$Rinv)
  delta.alpha = prior$delta.alpha; delta.beta = prior$delta.beta

  n0 <- length(y.Hist); n1 <- length(y.Cur)
  yVec0 = matrix(y.Hist, ncol = 1)
  yVec1 = matrix(y.Cur, ncol = 1)

  if(is.null(x.Cur) && is.null(x.Hist)){
    x.Cur = matrix(1, ncol = 1, nrow = n1)
    x.Hist = matrix(1, ncol = 1, nrow = n0)
  }

  if(is.data.frame(x.Cur)) x.Cur = as.matrix(x.Cur)
  if(is.data.frame(x.Hist)) x.Hist = as.matrix(x.Hist)

  x.Cur = cbind(1, x.Cur); x.Hist = cbind(1, x.Hist)

  if(nrow(x.Cur) != n1){stop("number of observations in covariates must match response")}
  if(nrow(x.Hist) != n0){stop("number of observations in covariates must match response")}
  if(ncol(x.Hist) != ncol(x.Cur)){stop("dimension of x.Hist must equal to x.Cur")}

  XtX0 = unname(t(x.Hist)%*%x.Hist)
  XtX1 = unname(t(x.Cur)%*%x.Cur)

  betaHat0 = solve(XtX0)%*%t(x.Hist)%*%yVec0
  betaHat1 = solve(XtX1)%*%t(x.Cur)%*%yVec1

  S0 = as.numeric(t(yVec0 - x.Hist%*%betaHat0)%*%(yVec0 - x.Hist%*%betaHat0))
  S1 = as.numeric(t(yVec1 - x.Cur%*%betaHat1)%*%(yVec1 - x.Cur%*%betaHat1))

  ncovariate = ncol(x.Cur)

  if(!(b == 0 | b == 1)){stop("b must be either 0 or 1")}
  if(b == 1){
    if(nrow(R) != ncovariate){stop("Dimension of R must equal to dimension of covariates")}
    if(nrow(mu0) != ncovariate){stop("Number of rows of mu0 must equal to dimension of covariates")}
  }
  #### Normalized Power Prior for Normal Log Marginal Posterior (unnormalized) of Delta
  xseq <- seq(.Machine$double.eps, 1-.Machine$double.eps, length = npoints)

  if(b == 1){
    LogPostNDelta <- function(x){
      betaS = solve(R+x*XtX0)%*%(R%*%mu0 + x*XtX0%*%betaHat0)
      H0 = as.numeric(t(mu0-betaHat0)%*%XtX0%*%(solve(R+x*XtX0))%*%R%*%(mu0-betaHat0))
      H1 = as.numeric(t(betaS-betaHat1)%*%XtX1%*%(solve(R+x*XtX0+XtX1))%*%(R+x*XtX0)%*%(betaS-betaHat1))
      lgM = (n1/2)*log(S1+H1+x*(S0+H0))+(a+x*n0/2-1)*log(1+(S1+H1)/(x*(S0+H0)))
      lden = (delta.alpha-1)*log(x)+(delta.beta-1)*log(1-x)+
        log(abs(det(R+x*XtX0)))/2+lgamma(a+n1/2+x*n0/2-1)-
        log(abs(det(R+x*XtX0+XtX1)))/2-lgamma(a+x*n0/2-1)-lgM
      return(lden)
    }
  }
  if(b == 0){
    deltamin = max(0, (ncovariate + 2 -2*a)/n0)
    LogPostNDelta <- function(x){
      H1 = as.numeric(t(betaHat0-betaHat1)%*%XtX1%*%(solve(x*XtX0+XtX1))%*%(x*XtX0)%*%(betaHat0-betaHat1))
      lgM = (n1/2)*log(S1+H1+x*(S0))+(a+x*n0/2-ncovariate/2-1)*log(1+(S1+H1)/(x*(S0)))
      lden = (delta.alpha-1)*log(x)+(delta.beta-1)*log(1-x)+
        log(abs(det(x*XtX0)))/2+lgamma(a+n1/2+x*n0/2-ncovariate/2-1)-
        log(abs(det(x*XtX0+XtX1)))/2-lgamma(a+x*n0/2-ncovariate/2-1)-lgM
      out = ifelse(x>deltamin, lden, log(.Machine$double.xmin))
      return(out)
    }
  }
  propDen = c()
  for(i in 1:length(xseq)){
  propDen[i] = LogPostNDelta(xseq[i])
  }
  NPPpostLik <- data.frame(x = xseq, logden = propDen)
  modedelta <- NPPpostLik[which.max(NPPpostLik$logden), 1]
  return(modedelta)
}
