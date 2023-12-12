LMMNPP_MCMC1 <- function(D0, X, Y, a0, b, mu0, R, delta_ini, prop_delta,
                         prior_delta_alpha, prior_delta_beta, prop_delta_alpha,
                         prop_delta_beta, rw_delta, nsample, burnin, thin){
  k = length(D0)
  l = ncol(X)
  n = length(Y)
  a = l/2+a0+1
  n0 = rep(0, k)
  X0X0 = list()
  X0Y0 = matrix(0, l, k)
  beta0hat = matrix(0, l, k)
  S0 = rep(0, k)
  beta0hatX0Y0 = rep(0, k)
  for(i in 1:k){
    n0[i] = nrow(D0[[i]])
    X0X0[[i]] = crossprod(D0[[i]][,1:l])
    X0Y0[, i] = crossprod(D0[[i]][,1:l], D0[[i]][,l+1])
    beta0hat[, i] = solve(X0X0[[i]])%*%X0Y0[, i]
    S0[i] = crossprod(D0[[i]][,l+1]-D0[[i]][,1:l]%*%beta0hat[, i])
    beta0hatX0Y0[i] = t(beta0hat[, i])%*%X0Y0[, i]
  }
  XtX = crossprod(X)
  XtY = crossprod(X, Y)
  betahat = solve(XtX)%*%crossprod(X, Y)
  S = crossprod(Y-X%*%betahat)
  #posterior distribution of delta
  logPostdelta = function(d){
    dX0X0 = array(0, dim=c(l, l, k))
    dX0Y0 = matrix(0, l, k)
    for(i in 1:k){
      dX0X0[, , i] = X0X0[[i]]*d[i]
      dX0Y0[, i] = X0Y0[, i]*d[i]
    }
    dX0X0_sum = apply(dX0X0, c(1, 2), sum, na.rm = TRUE)
    dX0Y0_sum = rowSums(dX0Y0)
    betastar = solve(R+dX0X0_sum)%*%(R%*%mu0+dX0Y0_sum)
    H0d = 2*b+sum(d*S0)+t(mu0)%*%R%*%mu0+sum(d*beta0hatX0Y0)-t(R%*%mu0+dX0Y0_sum)%*%betastar
    Hd = H0d+S+
         t(betastar-betahat)%*%XtX%*%solve(R+dX0X0_sum+XtX)%*%(R+dX0X0_sum)%*%(betastar-betahat)
    nu0 = (sum(d*n0)-l)/2+a-1
    out = 0.5*(log(det(R+dX0X0_sum))-log(det(R+dX0X0_sum+XtX)))+log(gamma(nu0+n/2))-
          log(gamma(nu0))+nu0*log(H0d/2)-(nu0+n/2)*log(Hd/2)+
          sum((prior_delta_alpha-1)*log(d)+(prior_delta_beta-1)*log(1-d))
    return(out)
  }
  if(is.null(delta_ini)){
    delta_ini = rep(1/k, k)
  }
  delta_sample = matrix(0, nrow = nsample, ncol = k)
  delta_cur = delta_sample[1, ] = delta_ini
  beta_sample = matrix(0, nrow = nsample, ncol = l)
  sigma_sample = rep(0, nsample)
  counter = 0
  niter = nsample*thin + burnin
  for (i in 2:niter){
    if(prop_delta == "IND"){
      delta_prop = rbeta(k, prop_delta_alpha, prop_delta_beta)
      llik.prop = logPostdelta(delta_prop)
      llik.cur = logPostdelta(delta_cur)
      logr = min(0, (llik.prop-llik.cur+
                     sum(dbeta(delta_cur,prop_delta_alpha,prop_delta_beta,log=TRUE)) -
                     sum(dbeta(delta_prop,prop_delta_alpha,prop_delta_beta,log=TRUE))))
    }
    if(prop_delta == "RW"){
      lgdelta_cur = log(delta_cur/(1-delta_cur))
      lgdelta_prop = rnorm(k, mean = lgdelta_cur, sd = rw_delta)
      delta_prop = exp(lgdelta_prop)/(1+exp(lgdelta_prop))
      llik.prop = logPostdelta(delta_prop)
      llik.cur = logPostdelta(delta_cur)
      logr = min(0, (llik.prop-llik.cur+
                     sum(log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur))))
    }
    if(runif(1) <= exp(logr)){
      delta_cur = delta_prop;
      counter = counter+1
    }
    if(i > burnin & (i-burnin) %% thin==0){
      delta_sample[(i-burnin)/thin, ] = delta_cur
    }
  }
  #Generate beta and sigma conditional on delta
  for(i in 1:nsample){
    deltaX0X0 = array(0, dim=c(l, l, k))
    deltaX0Y0 = matrix(0, l, k)
    for(j in 1:k){
      deltaX0X0[, , j] = X0X0[[j]]*delta_sample[i, j]
      deltaX0Y0[, j] = X0Y0[, j]*delta_sample[i, j]
    }
    deltaX0X0_sum = apply(deltaX0X0, c(1, 2), sum, na.rm = TRUE)
    deltaX0Y0_sum = rowSums(deltaX0Y0)
    betastar = solve(R+deltaX0X0_sum)%*%(R%*%mu0+deltaX0Y0_sum)
    betatilde = solve(R+deltaX0X0_sum+XtX)%*%(R%*%mu0+deltaX0Y0_sum+XtY)
    Hd = 2*b+sum(delta_sample[i, ]*S0)+t(mu0)%*%R%*%mu0+sum(delta_sample[i, ]*beta0hatX0Y0)-
         t(R%*%mu0+deltaX0Y0_sum)%*%betastar+S+
         t(betastar-betahat)%*%XtX%*%solve(R+deltaX0X0_sum+XtX)%*%(R+deltaX0X0_sum)%*%(betastar-betahat)
    nu0 = (sum(delta_sample[i, ]*n0)-l)/2+a-1
    beta_sample[i, ] = rmt(1, mean = betatilde, S = as.numeric(Hd/(2*nu0+n))*solve(R+deltaX0X0_sum+XtX), df = 2*nu0+n)
    sigma_sample[i] = r_igamma(1, shape = nu0+n/2, rate = Hd/2)
  }
  out = list(acceptrate = counter/niter, beta = beta_sample,
             sigma = sigma_sample, delta = delta_sample)
  return(out)
}
