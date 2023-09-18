PoiMNPP_MCMC2 <- function(n0,n,prior_lambda,prop_delta,
                          prior_delta_alpha,prior_delta_beta,rw_delta,
                          delta_ini,nsample,burnin,thin){
  k = length(n0)
  LogPostPoidelta = function(x){
    gamman0 = x*n0
    out = sum(lgamma(gamman0+n+prior_lambda[1])-lgamma(gamman0+prior_lambda[1])+
                (gamman0+prior_lambda[1])*log(x+prior_lambda[2])- 
                (gamman0+prior_lambda[1]+n)*log(x+prior_lambda[2]+1))+
      sum((prior_delta_alpha-1)*log(x)+(prior_delta_beta-1)*log(1-x))
    return(out)
  }
  if(is.null(delta_ini)){
    delta_ini = rep(1/k,k)
  }
  delta_cur = delta_ini
  delta_sample = matrix(0, nrow = nsample, ncol = k)
  counter = 0
  niter = nsample*thin + burnin
  
  for (i in 2:niter){
    if(prop_delta == "IND"){
      delta_prop = rbeta(k, prior_delta_alpha, prior_delta_beta)
      llik.prop = LogPostPoidelta(delta_prop)
      llik.cur = LogPostPoidelta(delta_cur)
      logr = min(0, (llik.prop-llik.cur+
                       sum(dbeta(delta_cur,prior_delta_alpha,prior_delta_beta,log=TRUE)) -
                       sum(dbeta(delta_prop,prior_delta_alpha,prior_delta_beta,log=TRUE))))
    }
    if(prop_delta == "RW"){
      lgdelta_cur = log(delta_cur/(1-delta_cur))
      lgdelta_prop = rnorm(k, mean = lgdelta_cur, sd = rw_delta)
      delta_prop = exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      llik.prop = LogPostPoidelta(delta_prop)
      llik.cur = LogPostPoidelta(delta_cur)
      logr = min(0, (llik.prop-llik.cur+
                       sum(log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur))))
    }
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop;
      counter = counter+1
    }
    if(i > burnin & (i-burnin) %% thin==0) {
      delta_sample[(i-burnin)/thin,] = delta_cur
    }
  }
  # Generate lambda conditional on delta 
  s1 = delta_sample%*%n0 + n + prior_lambda[1]
  s2_inv = rowSums(delta_sample) + 1 + prior_lambda[2]
  s2=1/s2_inv
  lambda = rgamma(nsample, s1, scale=s2)#scale is 1/s
  ans = list(acceptrate = counter/niter, lambda = lambda, delta = delta_sample)
  return(ans)
}