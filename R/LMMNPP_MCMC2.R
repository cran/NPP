LMMNPP_MCMC2 <- function(D0, X, Y, a0, b, mu0, R, delta_ini, prop_delta,
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
    lambda0 = array(0, dim=c(l, l, k))
    lambda0_det = rep(0, k)
    dX0Y0 = matrix(0, l, k)
    belambe = rep(0, k)
    H0d = rep(0, k)
    nu0 = rep(0, k)
    for(i in 1:k){
      lambda0[, , i] = R+X0X0[[i]]*d[i]
      lambda0_det[i] = det(lambda0[, , i])
      dX0Y0[, i] = X0Y0[, i]*d[i]
      belambe[i] = t(R%*%mu0+X0Y0[, i]*d[i])%*%solve(lambda0[, , i])%*%(R%*%mu0+X0Y0[, i]*d[i])
      H0d[i] = 2*b+d[i]*S0[i]+
               d[i]*t(mu0-beta0hat[, i])%*%X0X0[[i]]%*%solve(lambda0[, , i])%*%R%*%(mu0-beta0hat[, i])
      nu0[i] = a-1+(n0[i]*d[i]-l)/2
    }
    dX0Y0_sum = rowSums(dX0Y0)
    lambda0_sum = apply(lambda0, c(1, 2), sum, na.rm = TRUE)
    H0 = sum(H0d)+sum(belambe)-t(k*R%*%mu0+dX0Y0_sum)%*%solve(lambda0_sum)%*%(k*R%*%mu0+dX0Y0_sum)
    H = H0+S+
        t(solve(lambda0_sum)%*%(k*R%*%mu0+dX0Y0_sum)-betahat)%*%XtX%*%solve(lambda0_sum+XtX)%*%lambda0_sum%*%(solve(lambda0_sum)%*%(k*R%*%mu0+dX0Y0_sum)-betahat)
    nu = (sum(d*n0)+n-l)/2+k*a-1
    out = sum(0.5*log(lambda0_det)+nu0*log(H0d/2)-log(gamma(nu0)))-
          0.5*log(det(lambda0_sum+XtX))+log(gamma(nu))-nu*log(H/2)+
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
      delta_cur = delta_prop
      counter = counter+1
    }
    if(i > burnin & (i-burnin) %% thin==0){
      delta_sample[(i-burnin)/thin, ] = delta_cur
    }
  }
  #Generate beta and sigma conditional on delta
  for(i in 1:nsample){
    lambda0 = array(0, dim=c(l, l, k))
    lambda0_det = rep(0, k)
    dX0Y0 = matrix(0, l, k)
    belambe = rep(0, k)
    H0d = rep(0, k)
    nu0 = rep(0, k)
    for(j in 1:k){
      lambda0[, , j] = R+X0X0[[j]]*delta_sample[i, j]
      lambda0_det[j] = det(lambda0[, , j])
      dX0Y0[, j] = X0Y0[, j]*delta_sample[i, j]
      belambe[j] = t(R%*%mu0+X0Y0[, j]*delta_sample[i, j])%*%solve(lambda0[, , j])%*%(R%*%mu0+X0Y0[, j]*delta_sample[i, j])
      H0d[j] = 2*b+delta_sample[i, j]*S0[j]+
               delta_sample[i, j]*t(mu0-beta0hat[, j])%*%X0X0[[j]]%*%solve(lambda0[, , j])%*%R%*%(mu0-beta0hat[, j])
      nu0[j] = a-1+(n0[j]*delta_sample[i, j]-l)/2

    }
    dX0Y0_sum = rowSums(dX0Y0)
    lambda0_sum = apply(lambda0, c(1, 2), sum, na.rm = TRUE)
    H0 = sum(H0d)+sum(belambe)-t(k*R%*%mu0+dX0Y0_sum)%*%solve(lambda0_sum)%*%(k*R%*%mu0+dX0Y0_sum)
    H = H0+S+
        t(solve(lambda0_sum)%*%(k*R%*%mu0+dX0Y0_sum)-betahat)%*%XtX%*%solve(lambda0_sum+XtX)%*%lambda0_sum%*%(solve(lambda0_sum)%*%(k*R%*%mu0+dX0Y0_sum)-betahat)
    betatilde = solve(lambda0_sum+XtX)%*%(k*R%*%mu0+dX0Y0_sum+XtY)
    nu = (sum(delta_sample[i, ]*n0)+n-l)/2+k*a-1
    beta_sample[i, ] = rmt(1, mean = betatilde, S = as.numeric(H/(2*nu))*solve(lambda0_sum+XtX), df = 2*nu)
    sigma_sample[i] = r_igamma(1, shape = nu, rate = H/2)
  }
  out = list(acceptrate = counter/niter, beta = beta_sample,
             sigma = sigma_sample, delta = delta_sample)
  return(out)
}
