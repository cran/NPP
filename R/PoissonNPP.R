
PoissonNPP_MCMC <- function(Data.Cur, Data.Hist, 
                            CompStat = list(n0 = NULL, mean0 = NULL, n1 = NULL, mean1 = NULL), 
                            prior = list(lambda.shape = 1/2, lambda.scale = 100, 
                                         delta.alpha = 1, delta.beta = 1), 
                            MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                            ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                            control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  
  # CompStat = list(mean0 = 10, n0 = 20, mean1 = 11, n1 = 30) 
  # prior = list(lambda.shape = 1/2, lambda.scale = 100,
  #              delta.alpha = 1, delta.beta = 1)
  # ind.delta.alpha= 1; ind.delta.beta= 1;
  # nsample = 10000; control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 1)
  
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
  
#### Normalized Power Prior for Poisson Log Marginal Posterior (unnormalized) of Delta
LogPostPoisDelta <- function(x){
  (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma(sum1+x*sum0+prior$lambda.shape)-
    lgamma(x*sum0+prior$lambda.shape)+
    (x*sum0+prior$lambda.shape)*log(n0*x+prior.lambda.rate)-
    (x*sum0+sum1+prior$lambda.shape)*log(n0*x+n1+prior.lambda.rate)
}


if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
 delta_cur <- delta.ini
 delta <- rep(delta.ini, nsample)
 counter <- 0
 niter <- nsample*control.mcmc$thin + control.mcmc$burnin
 
 
 for (i in 1:niter){
  
   
  ### Update delta with RW MH for Logit delta
  if(MCMCmethod == 'RW'){
     lgdelta_cur <- log(delta_cur/(1-delta_cur))
     lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
     delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))

     llik.prop <- LogPostPoisDelta(delta_prop)
     llik.cur <- LogPostPoisDelta(delta_cur)
  
     logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
  }
   
  if(MCMCmethod == 'IND'){
     delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
     llik.prop <- LogPostPoisDelta(delta_prop)
     llik.cur <- LogPostPoisDelta(delta_cur)
     
     logr <- min(0, (llik.prop-llik.cur+
                     dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                     dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))
  }
   
  if(runif(1) <= exp(logr)){
    delta_cur = delta_prop; counter = counter+1
  }
   
  if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0) {
     delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
   }
}

# Generate lambda conditional on delta; vectorized
lambda <- rgamma(nsample, shape = sum1+delta*sum0+prior$lambda.shape, 
                          rate = n1+n0*delta+prior.lambda.rate)

# Calculate DIC
meanlambda <- mean(lambda)
D <- -2*(-n1*lambda + log(lambda)*sum1)
Dlambdabar <- -2*(-n1*meanlambda + log(meanlambda)*sum1)
DIC <- 2*mean(D)-Dlambdabar


return(list(lambda = lambda, delta = delta, acceptrate = counter/niter, DIC = DIC))
}
