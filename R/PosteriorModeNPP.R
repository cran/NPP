
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
                             prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1))
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
  propDen <-  (x*n0/2+prior$joint.a+prior$delta.alpha-2)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma((x*n0+n1-3)/2+prior$joint.a)-lgamma((x*n0-3)/2+prior$joint.a)-
    ((x*n0+n1-3)/2+prior$joint.a)*log((x*n1*(mean0-mean1)^2)/((x*n0+n1)*var0)+x+(n1*var1)/(n0*var0))
  
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


  
  