
#### Function to sample from location-scale t distribution
rt_ls <- function(n, df, location, scale) rt(n,df)*scale + location

#### Function to sample from Inverse Gamma
r_igamma <- function(n, shape, rate = 1, scale = 1/rate){
  if(missing(rate) && !missing(scale)) rate <- 1/scale
  1/rgamma(n, shape, rate)
}

#### Assume Joint Prior of mu and sigmasq is (1/sigmasq)^a
#### a = 1.5 Jeffreys; a = 1 reference
#### When provide CompStat, then ignore data, when not provide, calculate from data
NormalNPP_MCMC <- function(Data.Cur, Data.Hist, 
                           CompStat = list(n0 = NULL, mean0 = NULL, var0 = NULL, n1 = NULL, mean1 = NULL, var1 = NULL), 
                           prior = list(joint.a = 1.5, 
                                        delta.alpha = 1, delta.beta = 1), 
                           MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                           ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000, 
                           control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
  mean0 <- mean(Data.Hist)
  n0 <- length(Data.Hist)
  var0 <- var(Data.Hist)*(n0-1)/n0
  
  mean1 <- mean(Data.Cur)
  n1 <- length(Data.Cur)
  var1 <- var(Data.Cur)*(n1-1)/n1
  }else{
    n0 <- CompStat$n0
    mean0 <- CompStat$mean0
    var0 <- CompStat$var0
    
    n1 <- CompStat$n1
    mean1 <- CompStat$mean1
    var1 <- CompStat$var1
  }
  #### Prior pi(mu,sigmasq) propto (1/sigmasq)^a
  ## a = 1 for Reference Prior; a = 1.5 for Jeffrey's Prior
  #### Normalized Power Prior for Normal Log Marginal Posterior (unnormalized) of Delta

    LogPostNDelta <- function(x){
    
    (x*n0/2+prior$joint.a+prior$delta.alpha-2)*log(x)+(prior$delta.beta-1)*log(1-x)+
    lgamma((x*n0+n1-3)/2+prior$joint.a)-lgamma((x*n0-3)/2+prior$joint.a)-
    ((x*n0+n1-3)/2+prior$joint.a)*log((x*n1*(mean0-mean1)^2)/((x*n0+n1)*var0)+x+(n1*var1)/(n0*var0))-
    log(x*n0+n1)/2
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
      
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostNDelta(delta_prop)
      llik.cur <- LogPostNDelta(delta_cur)
      
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
  
#### vectorized generate mu and sigmasq conditional on delta 
  K = (delta*n0*n1*(mean0-mean1)^2/(delta*n0+n1) + delta*n0*var0 + n1*var1 )/2
  mu = rt_ls(nsample, df= delta*n0+n1+2*prior$joint.a-3, location= (delta*n0*mean0+n1*mean1)/(delta*n0 + n1), 
                 scale= sqrt(2*K/( (delta*n0+n1+2*prior$joint.a-3)*(delta*n0+n1) )) )
  
  sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1+2*prior$joint.a-3)/2, rate = K)
  
  meanmu <- mean(mu)
  meansigmasq <- mean(sigmasq)
  
  ### DIC without constant term
  D <- -2*(-n1*log(sigmasq)/2 -n1*(var1 + (mu-mean1)^2 )/(2*sigmasq))
  Dbar <- -2*(-n1*log(meansigmasq)/2 -n1*(var1 + (meanmu-mean1)^2 )/(2*meansigmasq))
  DIC <- 2*mean(D)-Dbar
  
  return(list(mu = mu, sigmasq = sigmasq, delta = delta, 
              acceptrate = counter/niter, DIC = DIC))
}


