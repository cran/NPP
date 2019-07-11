r_dir<-function(n,alpha)
  ## generate n random deviates from the Dirichlet function with shape
  ## parameters alpha is a vector of shapes
  ## source: gtools
{
  l<-length(alpha);
  x<-matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE);
  sm<-x%*%rep(1,l);
  x/as.vector(sm);
}

## Default: Jeffrey's Prior for theta
MultinomialNPP_MCMC <- function(Data.Cur = c(10, 10, 10), Data.Hist = c(10, 10, 10),  
                        CompStat = list(n0 = NULL, n1 = NULL), 
                        prior = list(theta.dir = c(0.5, 0.5, 0.5), 
                                     delta.alpha = 1, delta.beta = 1), 
                        MCMCmethod = 'IND', rw.logit.delta = 0.1, 
                        ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000, 
                        control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{
  if(missing(CompStat)){
    if(length(Data.Cur) != length(Data.Cur)) stop("Number of Categories Must Be Equal!")
    k <- length(Data.Cur)
    n0 <- Data.Hist
    n1 <- Data.Cur
  }else{
    if(length(CompStat$n0) != length(CompStat$n1)) stop("Number of Categories Must Be Equal!")
    k <- length(CompStat$n1)
    n0 <- CompStat$n0
    n1 <- CompStat$n1
  }
  
  n0sum <- sum(n0); n1sum <- sum(n1)
  alphasum <- sum(prior$theta.dir)
  
  #### Normalized Power Prior for Dirichlet Log Marginal Posterior (unnormalized) of Delta
  LogPostDirDelta <- function(x){
    s1 <- sum(lgamma(n0*x+n1+prior$theta.dir))
    s2 <- sum(lgamma(n0*x+prior$theta.dir))
    result <- ((prior$delta.alpha-1)*log(x)+(prior$delta.beta-1)*log(1-x)+
      lgamma(alphasum+x*n0sum)+s1-s2-
      lgamma(alphasum+x*n0sum+n1sum))
    return(result)
  }
  
  if(is.null(control.mcmc$delta.ini)) delta.ini = 0.5
  delta_cur <- delta.ini
  delta <- rep(delta.ini, nsample)
  theta <- matrix(0, ncol = k, nrow = nsample)
  counter <- 0
  niter <- nsample*control.mcmc$thin + control.mcmc$burnin
  
  for (i in 1:niter){
    ### Update delta with RW MH for Logit delta
    
    if(MCMCmethod == 'RW'){
      lgdelta_cur <- log(delta_cur/(1-delta_cur))
      lgdelta_prop <- rnorm(1, mean = lgdelta_cur, sd = sqrt(rw.logit.delta))
      delta_prop <- exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      
      llik.prop <- LogPostDirDelta(delta_prop)
      llik.cur <- LogPostDirDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
    }
    
    if(MCMCmethod == 'IND'){
      delta_prop <- rbeta(1, shape1 = ind.delta.alpha, shape2 = ind.delta.beta)
      llik.prop <- LogPostDirDelta(delta_prop)
      llik.cur <- LogPostDirDelta(delta_cur)
      
      logr <- min(0, (llik.prop-llik.cur+
                        dbeta(delta_cur, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE) -
                        dbeta(delta_prop, shape1 = ind.delta.alpha, shape2 = ind.delta.beta, log = TRUE)))
    }
    
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop; counter = counter+1
    }
    
    if( i > control.mcmc$burnin & (i-control.mcmc$burnin)%%control.mcmc$thin==0)
    {
      delta[(i-control.mcmc$burnin)/control.mcmc$thin] <- delta_cur
      theta[(i-control.mcmc$burnin)/control.mcmc$thin, ] <- r_dir(1, alpha = n0*delta_cur+n1+prior$theta.dir)  
    }
  }
  
  # Calculate DIC
    meantheta <- apply(theta, 2, mean)
    D <- -2*( log(theta)%*%as.vector(n1))
    Dpbar <- -2*(sum(n1*log(meantheta)))
    DIC <- 2*mean(D)-Dpbar
  
  return(list(theta = theta, delta = delta, acceptrate = counter/niter, DIC = DIC))
}
