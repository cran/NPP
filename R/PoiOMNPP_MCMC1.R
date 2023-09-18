PoiOMNPP_MCMC1 <- function(n0,n,prior_gamma,prior_lambda,
                           gamma_ind_prop,gamma_ini,nsample,burnin,thin){
  node = c(nsample/5, nsample/2)
  k = length(n0)
  LogPostPoigamma = function(x){
    cumsumx = cumsum(x[1:k])
    gamman0 = sum(cumsumx*n0)
    gamma0 = sum(cumsumx)
    out = lgamma(gamman0+n+prior_lambda[1])-lgamma(gamman0+prior_lambda[1])+
      (gamman0+prior_lambda[1])*log(gamma0+prior_lambda[2]) -
      (gamman0+prior_lambda[1]+n)*log(gamma0+prior_lambda[2]+1)+sum((prior_gamma-1)*log(x))
    return(out)
  }
  if(is.null(gamma_ini)){
    gamma_ini = rep(1/(k+1),k+1)
  }
  gamma_cur = gamma_ini
  gamma_sample = matrix(0, nrow = nsample, ncol = k+1)
  counter = 0
  niter = nsample*thin + burnin

  for (i in 1:niter){
    ## Independent Dirichlet Proposal with adjustment at selected iterations
    ii = (i-burnin)/thin
    if(ii %in% node){
      gamma_ind_prop = c(dirichmomcopy(gamma_sample[1:ii, ]))
    }
    gamma_prop = rdirichletcopy(1, alpha = gamma_ind_prop)
    lpost_prop = LogPostPoigamma(gamma_prop)
    lpost_cur = LogPostPoigamma(gamma_cur)
    logr = min(0, (lpost_prop-lpost_cur+
                     log(ddirichletcopy(gamma_cur, alpha = gamma_ind_prop)) -
                     log(ddirichletcopy(gamma_prop, alpha = gamma_ind_prop))))

    if(runif(1) <= exp(logr)){
      gamma_cur = gamma_prop;
      counter = counter+1
    }
    if( i > burnin & (i-burnin) %% thin==0) {
      gamma_sample[ii, ] = gamma_cur
    }
  }
  gamma_minus = gamma_sample[, 1:k]
  delta_sample = t(apply(gamma_minus, 1, cumsum))
  # Generate lambda conditional on gamma; vectorized
  s1 = delta_sample%*%n0 + n + prior_lambda[1]
  s2_inv = rowSums(delta_sample) + 1 + prior_lambda[2]
  s2=1/s2_inv
  lambda = rgamma(nsample, s1, scale=s2)#scale is 1/s
  ans = list(acceptrate = counter/niter, lambda = lambda, delta = delta_sample)
  return(ans)
}
