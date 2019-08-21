#### Function to sample from location-scale t and multivariate t distribution
rt_ls <- function(n, df, location, scale) rt(n,df)*scale + location

rmnorm <- function(n=1, mean, sqrt=NULL)
{
  sqrt.varcov <- sqrt
  d <- if(is.matrix(sqrt.varcov)) ncol(sqrt.varcov) else 1
  mean <- outer(rep(1,n), as.vector(matrix(mean,d)))
  drop(mean + t(matrix(rnorm(n*d), d, n)) %*% sqrt.varcov)
}

rmt <- function(n=1, mean, S, df=Inf, sqrt=NULL)
{
  sqrt.S <- if(is.null(sqrt)) chol(S) else sqrt
  d <- if(is.matrix(sqrt.S)) ncol(sqrt.S) else 1
  x <- if(df==Inf) 1 else rchisq(n, df)/df
  z <- rmnorm(n, mean = rep(0, d), sqrt=sqrt.S)
  mean <- outer(rep(1, n), as.vector(matrix(mean,d)))
  drop(mean + z/sqrt(x))
}

#### Function to sample from Inverse Gamma
#### The rate here is rate for Gamma Distribution, in Inverse Gamma it is scale
r_igamma <- function(n, shape, rate = 1, scale = 1/rate){
  if(missing(rate) && !missing(scale)) rate <- 1/scale
  1/rgamma(n, shape, rate)
}

#### Assume Joint Prior of beta and sigmasq is (1/sigmasq)^a*Normal(mu0, sigmasq*R^-1) when b = 1
#### Assume Joint Prior of beta and sigmasq is (1/sigmasq)^a when b = 0
LMNPP_MCMC <- function(y.Cur, y.Hist, x.Cur = NULL, x.Hist = NULL,
                       prior = list(a = 1.5, b = 0, mu0 = 0,
                                    Rinv = matrix(1, nrow = 1),
                                    delta.alpha = 1, delta.beta = 1),
                       MCMCmethod = 'IND', rw.logit.delta = 0.1,
                       ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000,
                       control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
{

  # prior = list(a = 1, b = 0, mu0 = matrix(0, nrow = 1),
  #              R = matrix(1, nrow = 1),
  #              delta.alpha = 1, delta.beta = 1)

  ### Start
  a = prior$a; b = prior$b
  mu0 = as.matrix(prior$mu0); R = solve(prior$Rinv)
  delta.alpha = prior$delta.alpha; delta.beta = prior$delta.beta

  n0 <- length(y.Hist); n1 <- length(y.Cur)
  yVec0 = matrix(y.Hist, ncol = 1)
  yVec1 = matrix(y.Cur, ncol = 1)

  if(is.null(x.Cur) && is.null(x.Hist)){
    x.Cur = matrix(1, ncol = 1, nrow = n1)
    x.Hist = matrix(1, ncol = 1, nrow = n0)
  }

  if(is.data.frame(x.Cur)) x.Cur = as.matrix(x.Cur)
  if(is.data.frame(x.Hist)) x.Hist = as.matrix(x.Hist)

  x.Cur = cbind(1, x.Cur); x.Hist = cbind(1, x.Hist)

  if(nrow(x.Cur) != n1){stop("number of observations in covariates must match response")}
  if(nrow(x.Hist) != n0){stop("number of observations in covariates must match response")}
  if(ncol(x.Hist) != ncol(x.Cur)){stop("dimension of x.Hist must equal to x.Cur")}

  XtX0 = unname(t(x.Hist)%*%x.Hist)
  XtX1 = unname(t(x.Cur)%*%x.Cur)

  betaHat0 = solve(XtX0)%*%t(x.Hist)%*%yVec0
  betaHat1 = solve(XtX1)%*%t(x.Cur)%*%yVec1

  S0 = as.numeric(t(yVec0 - x.Hist%*%betaHat0)%*%(yVec0 - x.Hist%*%betaHat0))
  S1 = as.numeric(t(yVec1 - x.Cur%*%betaHat1)%*%(yVec1 - x.Cur%*%betaHat1))

  ncovariate = ncol(x.Cur)

  if(!(b == 0 | b == 1)){stop("b must be either 0 or 1")}
  if(b == 1){
    if(nrow(R) != ncovariate){stop("Dimension of R must equal to dimension of covariates")}
    if(nrow(mu0) != ncovariate){stop("Number of rows of mu0 must equal to dimension of covariates")}
  }
  #### Normalized Power Prior for Normal Log Marginal Posterior (unnormalized) of Delta

  if(b == 1){
    LogPostNDelta <- function(x){
      betaS = solve(R+x*XtX0)%*%(R%*%mu0 + x*XtX0%*%betaHat0)
      H0 = as.numeric(t(mu0-betaHat0)%*%XtX0%*%(solve(R+x*XtX0))%*%R%*%(mu0-betaHat0))
      H1 = as.numeric(t(betaS-betaHat1)%*%XtX1%*%(solve(R+x*XtX0+XtX1))%*%(R+x*XtX0)%*%(betaS-betaHat1))
      lgM = (n1/2)*log(S1+H1+x*(S0+H0))+(a+x*n0/2-1)*log(1+(S1+H1)/(x*(S0+H0)))
      lden = (delta.alpha-1)*log(x)+(delta.beta-1)*log(1-x)+
        log(abs(det(R+x*XtX0)))/2+lgamma(a+n1/2+x*n0/2-1)-
        log(abs(det(R+x*XtX0+XtX1)))/2-lgamma(a+x*n0/2-1)-lgM
      return(lden)
    }
  }
  if(b == 0){
    deltamin = max(0, (ncovariate + 2 -2*a)/n0)
    LogPostNDelta <- function(x){
      H1 = as.numeric(t(betaHat0-betaHat1)%*%XtX1%*%(solve(x*XtX0+XtX1))%*%(x*XtX0)%*%(betaHat0-betaHat1))
      lgM = (n1/2)*log(S1+H1+x*(S0))+(a+x*n0/2-ncovariate/2-1)*log(1+(S1+H1)/(x*(S0)))
      lden = (delta.alpha-1)*log(x)+(delta.beta-1)*log(1-x)+
        log(abs(det(x*XtX0)))/2+lgamma(a+n1/2+x*n0/2-ncovariate/2-1)-
        log(abs(det(x*XtX0+XtX1)))/2-lgamma(a+x*n0/2-ncovariate/2-1)-lgM
      out = ifelse(x>deltamin, lden, log(.Machine$double.xmin))
      return(out)
    }
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

  if(ncovariate == 1){
    mean0 <- mean(y.Hist); mean1 <- mean(y.Cur)
    var0 <- var(y.Hist)*(n0-1)/n0
    var1 <- var(y.Cur)*(n1-1)/n1
    K = (delta*n0*n1*(mean0-mean1)^2/(delta*n0+n1) + delta*n0*var0 + n1*var1)/2
    mu = rt_ls(nsample, df= delta*n0+n1+2*a-3, location= (delta*n0*mean0+n1*mean1)/(delta*n0 + n1),
               scale= sqrt(2*K/((delta*n0+n1+2*a-3)*(delta*n0+n1))))

    sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1+2*a-3)/2, rate = K)
    meanmu <- mean(mu); meansigmasq <- mean(sigmasq)

    ### DIC without constant term
    D <- -2*(-n1*log(sigmasq)/2 -n1*(var1 + (mu-mean1)^2)/(2*sigmasq))
    Dbar <- -2*(-n1*log(meansigmasq)/2 -n1*(var1 + (meanmu-mean1)^2)/(2*meansigmasq))
    DIC <- 2*mean(D)-Dbar

    return(list(beta = mu, sigmasq = sigmasq, delta = delta, acceptrate = counter/niter, DIC = DIC))
  }

  #### Sample the beta (non-vectorized input due to dimensionality) and then sigmasq (vectorized)
  if(ncovariate > 1){
    scalevec = rep(1, nsample)
    beta = matrix(1, nrow = nsample, ncol = ncovariate)

    if(b == 1){
      for(i in 1:nsample){
        betaS = solve(R+delta[i]*XtX0)%*%(R%*%mu0 + delta[i]*XtX0%*%betaHat0)
        mu = solve(R+delta[i]*XtX0+XtX1)%*%((R+delta[i]*XtX0)%*%betaS +XtX1%*%betaHat1)
        dfv = delta[i]*n0+n1+2*a - 2

        H0 = as.numeric(t(mu0-betaHat0)%*%XtX0%*%(solve(R+delta[i]*XtX0))%*%R%*%(mu0-betaHat0))
        H1 = as.numeric(t(betaS-betaHat1)%*%XtX1%*%(solve(R+delta[i]*XtX0+XtX1))%*%
                          (R+delta[i]*XtX0)%*%(betaS-betaHat1))

        Sigma = ((H1+S1+delta[i]*(S0+H0))/dfv)*(solve(R+delta[i]*XtX0 + XtX1))
        beta[i, ] = rmt(n = 1, mean = mu, S = Sigma, df = dfv)
        scalevec[i] = (H1+S1+delta[i]*(S0+H0))/2
      }
      sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1)/2+a-1, rate = scalevec)
    }

    if(b == 0){
      for(i in 1:nsample){
        mu = solve(delta[i]*XtX0 +XtX1)%*%(delta[i]*XtX0%*%betaHat0 +XtX1%*%betaHat1)
        dfv = delta[i]*n0+n1+2*a-ncovariate-2
        H1 = as.numeric(t(betaHat0-betaHat1)%*%XtX1%*%(solve(delta[i]*XtX0+XtX1))%*%
                          (delta[i]*XtX0)%*%(betaHat0-betaHat1))
        Sigma = ((H1+S1+delta[i]*S0)/dfv)*(solve(delta[i]*XtX0 + XtX1))
        beta[i, ] = rmt(n = 1, mean = mu, S = Sigma, df = dfv)
        scalevec[i] = (H1+S1+delta[i]*S0)/2
      }
      sigmasq= r_igamma(n = nsample, shape = (delta*n0+n1-ncovariate)/2+a-1, rate = scalevec)
    }
    yHat = beta%*%t(x.Cur)
    meansigmasq <- mean(sigmasq)
    yHatBar = apply(yHat, 2, mean)

    ### DIC without constant term

    SQdiffSum = apply((t(as.numeric(yVec1) - t(yHat)))^2, 1, sum)

    D <- -2*( -n1*log(sigmasq)/2 -SQdiffSum/(2*sigmasq) )
    Dbar <- -2*( -n1*log(meansigmasq)/2 -(sum((y.Cur-yHatBar)^2))/(2*meansigmasq) )
    DIC <- 2*mean(D)-Dbar

    return(list(beta = beta, sigmasq = sigmasq, delta = delta,
                acceptrate = counter/niter, DIC = DIC))
  }
}
