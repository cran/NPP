
BerNPP_MCMC <- function(Data.Cur = c(100, 50), Data.Hist = c(100, 50), 
                        CompStat = list(n0 = NULL, y0 = NULL, n1 = NULL, y1 = NULL), 
                         prior = list(p.alpha = 1, p.beta = 1, 
                                     delta.alpha = 1, delta.beta = 1), 
                        MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                        ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000, 
                        control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
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
  #### Normalized Power Prior for Bernoulli Log Marginal Posterior (unnormalized) of Delta
  
  LogPostBerDelta <- function(x){
    lbeta(x*y0 + y1 + prior$p.alpha, x*(n0-y0)+n1-y1+prior$p.beta) -
    lbeta(x*y0 + prior$p.alpha, x*(n0-y0) + prior$p.beta) + 
    (prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)
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
      
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostBerDelta(delta_prop)
      llik.cur <- LogPostBerDelta(delta_cur)
      
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
  
  # Generate p conditional on delta; vectorized 
  p <- rbeta(nsample, shape1 = delta*y0+y1+prior$p.alpha, 
                      shape2 = delta*(n0-y0)+(n1-y1)+prior$p.beta)    
  
  # Calculate DIC
  meanp <- mean(p)
  D <- -2*(y1*log(p)+(n1-y1)*log(1-p))
  Dpbar <- -2*(y1*log(meanp)+(n1-y1)*log(1-meanp))
  DIC <- 2*mean(D)-Dpbar
  ans <- list(p = p, delta = delta, acceptrate = counter/niter, DIC = DIC)
  
  class(ans) <- "NPP"
  return(ans)
}

