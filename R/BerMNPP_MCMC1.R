BerMNPP_MCMC1 <- function(n0,y0,n,y,prior_p,prior_delta_alpha,prior_delta_beta,prop_delta_alpha,
                          prop_delta_beta,delta_ini,prop_delta,rw_delta,nsample,burnin,thin){
  k = length(n0) #number of historical datasets
  LogPostBerdelta = function(x){
    deltay0 = sum(x*y0)
    deltanmy0 = sum(x*(n0-y0))
    out = lbeta(deltay0+y+prior_p[1],n-y+deltanmy0+prior_p[2])- 
      lbeta(deltay0+prior_p[1],deltanmy0+prior_p[2])+ 
      sum((prior_delta_alpha-1)*log(x)+(prior_delta_beta-1)*log(1-x))
    return(out)
  }
  if(is.null(delta_ini)){ 
    delta_ini = rep(1/k,k)
  }
  delta_sample = matrix(0, nrow = nsample, ncol = k)
  delta_cur = delta_sample[1,] = delta_ini
  counter = 0
  niter = nsample*thin+burnin
  
  for (i in 2:niter){
    if(prop_delta == "IND"){
      delta_prop = rbeta(k, prop_delta_alpha, prop_delta_beta)
      llik.prop = LogPostBerdelta(delta_prop)
      llik.cur = LogPostBerdelta(delta_cur)
      logr = min(0, (llik.prop-llik.cur+
                       sum(dbeta(delta_cur,prop_delta_alpha,prop_delta_beta,log=TRUE)) -
                       sum(dbeta(delta_prop,prop_delta_alpha,prop_delta_beta,log=TRUE))))
    }
    if(prop_delta == "RW"){
      lgdelta_cur = log(delta_cur/(1-delta_cur))
      lgdelta_prop = rnorm(k, mean = lgdelta_cur, sd = rw_delta)
      delta_prop = exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      llik.prop = LogPostBerdelta(delta_prop)
      llik.cur = LogPostBerdelta(delta_cur)
      logr = min(0, (llik.prop-llik.cur+
                       sum(log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur))))
    }
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop;
      counter = counter+1
    }
    if(i > burnin & (i-burnin) %% thin==0){
      delta_sample[(i-burnin)/thin,] = delta_cur
    }
  }
  # Generate p conditional on gamma; vectorized 
  s1 = delta_sample%*%y0 + y + prior_p[1]
  s2 = delta_sample%*%(n0-y0) + n - y + prior_p[2]
  p = rbeta(nsample, s1, s2)
  ans = list(acceptrate=counter/niter,p=p,delta=delta_sample)
  return(ans)
}