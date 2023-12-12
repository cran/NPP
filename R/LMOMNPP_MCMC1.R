LMOMNPP_MCMC1 <- function(D0, X, Y, a0, b, mu0, R, gamma_ini, prior_gamma,
                          gamma_ind_prop, nsample, burnin, thin, adjust){
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
  logPostgamma = function(x){
    d = cumsum(x[1:k])
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
          log(gamma(nu0))+nu0*log(H0d/2)-(nu0+n/2)*log(Hd/2)+sum((prior_gamma-1)*log(x))
    return(out)
  }
  if(is.null(gamma_ini)){
    gamma_ini = rep(1/(k+1), k+1)
  }
  gamma_cur = gamma_ini
  gamma_sample = matrix(0, nrow = nsample, ncol = k+1)
  beta_sample = matrix(0, nrow = nsample, ncol = l)
  sigma_sample = rep(0, nsample)
  counter = 0
  niter = nsample*thin + burnin
  node = c(nsample/4, nsample/2)
  for (i in 1:niter){
    ## Independent Dirichlet Proposal with adjustment at selected iterations
    ii = (i-burnin)/thin
    if(isTRUE(adjust)){
      if(ii %in% node){
        gamma_ind_prop = c(dirichmomcopy(gamma_sample[1:ii, ]))
      }
    }
    gamma_prop = rdirichletcopy(1, alpha = gamma_ind_prop)
    lpost_prop = min(logPostgamma(gamma_prop), 1e100)
    lpost_prop = max(lpost_prop, -1e100)                           ## Lazy way To Avoid Overflow
    lpost_cur = min(logPostgamma(gamma_cur), 1e100)
    lpost_cur = max(lpost_cur, -1e100)
    safelprior_cur = min(log(ddirichletcopy(gamma_cur, alpha = gamma_ind_prop)), 1e100)
    safelprior_cur = max(safelprior_cur, -1e100)
    safelprior_prop = min(log(ddirichletcopy(gamma_prop, alpha = gamma_ind_prop)), 1e100)
    safelprior_prop = max(safelprior_prop, -1e100)
    logr = min(0, (lpost_prop-lpost_cur+safelprior_cur-safelprior_prop))

    if(runif(1) <= exp(logr)){
      gamma_cur = gamma_prop
      counter = counter+1
    }
    if(i > burnin & (i-burnin) %% thin==0) {
      gamma_sample[ii, ] = gamma_cur
    }
  }
  gamma_minus = gamma_sample[, 1:k]
  delta_sample = t(apply(gamma_minus, 1, cumsum))
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
