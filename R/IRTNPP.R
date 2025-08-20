makePositiveDefinite = function(m, tol) {
  if (!is.matrix(m)) {
    m = as.matrix(m)
  }
  d = dim(m)[1]
  if ( dim(m)[2] != d ) {
    stop("Input matrix is not square!")
  }
  es = eigen(m)
  esv = es$values
  if (missing(tol)) {
    tol = d*max(abs(esv))*.Machine$double.eps
  }
  delta =  2*tol
  tau = pmax(0, delta - esv)
  dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
  return( m + dm )
}

Metro_Hastings = function(li_func,pars,prop_sigma=NULL,par_names=NULL, 
                          iterations=50000,burn_in=1000,adapt_par=c(100,20,0.5,0.75),quiet=TRUE,...){
  if (!is.finite(li_func(pars, ...))) 
    stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")
  if (is.null(par_names))
    par_names<-letters[1:length(pars)]
  if(!is.null(dim(prop_sigma)))
  {
    if( ( dim(prop_sigma)[1] != length(pars) ||  dim(prop_sigma)[2] != length(pars) ) && !is.null(prop_sigma) )
      stop("prop_sigma not of dimension length(pars) x length(pars)")
  }
  if(is.null(prop_sigma)) 
  {
    if(length(pars)!=1)
    {
      fit<-optim(pars,li_func,control=list("fnscale"=-1),hessian=TRUE,...)
      fisher_info<-solve(-fit$hessian)
      prop_sigma<-sqrt(diag(fisher_info))
      prop_sigma<-diag(prop_sigma)
    }else{
      prop_sigma<-1+pars/2
    }
  }
  prop_sigma<-makePositiveDefinite(prop_sigma)
  mu<-pars
  pi_X<-li_func(pars,...)           
  k_X<-pars
  trace<-array(dim=c(iterations,length(pars)))
  deviance<-array(dim=iterations) 
  announce<-floor(seq(iterations/10,iterations,length.out=10))
  for(i in 1:iterations)
  {
    k_Y<- MASS::mvrnorm(1,mu=k_X,Sigma=prop_sigma)
    pi_Y<-li_func(k_Y,...)    
    a_X_Y = (pi_Y)-(pi_X)          
    if(is.nan(a_X_Y))
      a_X_Y<--Inf                          
    if( log(runif(1,0,1)) <= a_X_Y)         
    {
      k_X = k_Y
      pi_X = pi_Y
    } 
    trace[i,]<-k_X 
    deviance[i]<-(-2*pi_X)                      
    if(i > adapt_par[1] && i %% adapt_par[2] == 0 && i < (adapt_par[4]*iterations) )       
    {   
      len<-floor(i*adapt_par[3]):i
      x<-trace[len,]
      N<-length(len)
      p_sigma <- (N-1) * var(x)/N
      p_sigma <-makePositiveDefinite(p_sigma)   
      if(!(0 %in% p_sigma) ) 
        prop_sigma<-p_sigma
    }
    if(!quiet && i %in% announce)
      print(paste('updating: ',i/iterations*100,'%',sep=''))
  }
  trace<-trace[burn_in:iterations,]
  DIC<-NULL
  D_bar<-mean(deviance[burn_in:iterations])
  if(length(pars)>1)
  {
    theta_bar<-sapply(1:length(pars),function(x){mean( trace[,x] )})
  }else
    theta_bar<-mean( trace )
  D_hat<-li_func(theta_bar,...)
  pD<-D_bar-D_hat
  DIC<-D_hat + 2*pD
  if(length(pars)>1)
    accept_rate<-length(unique(trace[,1]))/(iterations-burn_in) else
      accept_rate<-length(unique(trace))/(iterations-burn_in)
  val<-list("trace"=trace,"prop_sigma"=prop_sigma,"par_names"=par_names,"DIC"=DIC,'acceptance_rate'=accept_rate)
  class(val)<-"MHposterior"
  return(val)
}


IRTNPP=function(y,dseq,prior_mu,prior_sd,MCsize,disa,difa1,difa2,
                cut,prior_beta,prior_delta,disb,difb1,difb2,
                prop_delta,rw_delta,rw_n_beta,rw_u_beta,ind_delta,
                prop_beta,n_sample,burnin,thin){
  Y0 <- y[1:nrow(disa)]
  y1 <- rep(0,length(Y0))
  y2 <- rep(0,length(Y0))
  for(i in 1:length(Y0)){
    if(Y0[i]==2){
      y1[i]=1
    }else{
      y1[i]=0
    }
    if(Y0[i]==3){
      y2[i]=1
    }else{
      y2[i]=0
    }
  }
  Y <- y[(nrow(disa)+1):(nrow(disa)+nrow(disb))]
  y3 <- rep(0,length(Y))
  y4 <- rep(0,length(Y))
  for(i in 1:length(Y)){
    if(Y[i]==2){
      y3[i]=1
    }else{
      y3[i]=0
    }
    if(Y[i]==3){
      y4[i]=1
    }else{
      y4[i]=0
    }
  }
  k=ncol(disb)
  set.seed(111)
  lgr_is_fun = function(dseq){   
    prior_var = prior_sd^2
    mu0 = rep(prior_mu,k)  
    sample_d0 = t(mvtnorm::rmvnorm(n=MCsize, mean=mu0, sigma=diag(prior_var, k)))
    li_reg = function(pars){
      lprior = sum(dnorm(pars,prior_mu,prior_sd,log=TRUE))
      log_L = sum((1-y1-y2)*log(1-exp(disa%*%pars+difa1)/(1+exp(disa%*%pars+difa1)))+
                    y1*log(exp(disa%*%pars+difa1)/(1+exp(disa%*%pars+difa1))-exp(disa%*%pars+difa2)/(1+exp(disa%*%pars+difa2)))+
                    y2*log(exp(disa%*%pars+difa2)/(1+exp(disa%*%pars+difa2))))
      return(lprior+log_L)
    }
    mcmc_r = Metro_Hastings(li_func=li_reg,pars=rep(0.1, k),iterations=MCsize+999,quiet=TRUE) 
    sample_d1 = t(do.call(cbind,mcmc_r[1]))
    meanbeta = apply(sample_d1,1,mean)
    llikFunMat = function(betaMat){
      tmat1 = exp(disa %*% betaMat + difa1)
      tmat2 = exp(disa %*% betaMat + difa2)
      tol11 = .Machine$double.xmin; tol22 = .Machine$double.xmax              
      tmat1[which(tmat1 >= tol22)] = tol22
      tmat1[which(tmat1 <= tol11)] = tol11
      tmat2[which(tmat2 >= tol22)] = tol22
      tmat2[which(tmat2 <= tol11)] = tol11
      pmat1 = tmat1/(1+tmat1)
      pmat2 = tmat2/(1+tmat2)
      pmat3 = pmat1-pmat2
      tol1 = 1e-7; tol2 = .Machine$double.xmin              
      pmat1[which(pmat1 > (1-tol1))] = 1-tol1
      pmat1[which(pmat1 < tol2)] = tol2
      pmat2[which(pmat2 > (1-tol1))] = 1-tol1
      pmat2[which(pmat2 < tol2)] = tol2
      pmat3[which(pmat3 > (1-tol1))] = 1-tol1
      pmat3[which(pmat3 < tol2)] = tol2
      loglikmat = t(1-y1-y2) %*% log(1-pmat1)+t(y1) %*% log(pmat3)+t(y2) %*% log(pmat2)
      return(loglikmat)
    }
    lpriorFunMat = function(betaMat){
      t = colSums((betaMat-mu0)^2)
      return(-log(prior_sd)-t/(2*prior_sd^2))
    }
    lnp = llikFunMat(sample_d1)
    lnprior = lpriorFunMat(sample_d1)
    integrand_fun = function(d){
      if(d <= cut){
        lnpb1 = llikFunMat(sample_d0)
        lnpb = d*lnpb1
        num = exp(lnpb-max(lnpb))
        denom = sum(num)
        W = num/denom 
      }else{
        sample_d = sample_d1/sqrt(d)+meanbeta*(1-1/sqrt(d))
        lnpriorb = lpriorFunMat(sample_d)
        lnpb1 = llikFunMat(sample_d)
        lnpb = d*lnpb1
        num = exp(lnpriorb+lnpb-lnprior-lnp-max(lnpriorb+lnpb-lnprior-lnp))
        denom = sum(num)
        W = num/denom
      }
      est = sum(lnpb1*W)
      return(est)
    }
    hd = lapply(dseq,integrand_fun)
    est = c()
    m = length(dseq)
    for(i in 1:m){
      est[i] = hd[[i]]
    }
    return(est)
  }
  
  Ctmp = lgr_is_fun(dseq=dseq)
  
  IS_logC = function(grid, logintgrand){
    n_g = length(grid)
    if(grid[1] == 0){
      out = c(0, cumsum(diff(grid)*(logintgrand[1:(n_g-1)]+logintgrand[2:n_g])/2))
    }else{
      r = (logintgrand[2]-logintgrand[1])/(grid[2]-grid[1])
      logintgrand0 = logintgrand[1]-r*grid[1]
      grid = c(0, grid)
      logintgrand = c(logintgrand0, logintgrand)
      out = cumsum(diff(grid)*(logintgrand[1:n_g]+logintgrand[2:(n_g+1)])/2)
    }
    return(cbind(grid, out))
  }
  
  logC_grid = IS_logC(grid=dseq, logintgrand=Ctmp)
  delta=logC_grid[,1]
  logC=logC_grid[,2]
  lr_NPP_mcmc <- function(logC_seq=list(delta, logC)){
    lprior0_beta = function(beta) -sum((beta-prior_beta$mu)^2)/(2*prior_beta$sd*prior_beta$sd)
    lprior0_delta = function(d) (prior_delta$alpha-1)*log(d) + (prior_delta$beta-1)*log(1-d)
    tol1 = 1e-7; tol2 = .Machine$double.xmin
    logCFun = function(d){                                     
      id_l = sum(d-logC_seq$delta > 0); id_u = id_l+1
      slope = (logC_seq$logC[id_u]-logC_seq$logC[id_l])/(logC_seq$delta[id_u]-logC_seq$delta[id_l])
      dif = d-logC_seq$delta[id_l]
      return(slope*dif + logC_seq$logC[id_l])
    }
    loglik0 = function(beta, d){
      tmat1 = exp(disa%*%beta + difa1)
      tmat2 = exp(disa%*%beta + difa2)
      p1 = tmat1/(1+tmat1)
      p2 = tmat2/(1+tmat2)
      p3 = p1-p2
      p1[which(p1 > (1-tol1))] = 1-tol1
      p1[which(p1 < tol2)] = tol2
      p2[which(p2 > (1-tol1))] = 1-tol1
      p2[which(p2 < tol2)] = tol2
      p3[which(p3 > (1-tol1))] = 1-tol1
      p3[which(p3 < tol2)] = tol2
      out = sum((1-y1-y2)*log(1-p1)+y1*log(p3)+y2*log(p2))*d
      return(out)
    }
    loglik = function(beta){
      tmat3 = exp(disb%*%beta + difb1)
      tmat4 = exp(disb%*%beta + difb2)
      p4 = tmat3/(1+tmat3)
      p5 = tmat4/(1+tmat4)
      p6 = p4-p5
      p4[which(p4 > (1-tol1))] = 1-tol1
      p4[which(p4 < tol2)] = tol2
      p5[which(p5 > (1-tol1))] = 1-tol1
      p5[which(p5 < tol2)] = tol2
      p6[which(p6 > (1-tol1))] = 1-tol1
      p6[which(p6 < tol2)] = tol2
      out = sum((1-y3-y4)*log(1-p4)+y3*log(p6)+y4*log(p5))
      return(out)
    }
    logPostDelta = function(d){                                
      lprior0_delta(d)+loglik0(beta_cur, d)-logCFun(d)
    }
    logPostBeta = function(beta){                               
      lprior0_beta(beta)+loglik0(beta, delta_cur)+loglik(beta)
    }
    # Matrix storing samples of the \beta parameter
    beta_chain = matrix(0, nrow=n_sample, ncol=k)
    delta_chain = rep(0, n_sample)
    beta_cur = beta_chain[1,] = rep(0.1, k)
    delta_cur = delta_chain[1] = 0.5
    counter1 = 0; counter2 = 0
    n_mcmc = n_sample*thin + burnin
    # Start MH algorithm
    for (i in 2:n_mcmc) {
      # Update delta
      if(prop_delta == "IND"){
        delta_prop = rbeta(1, ind_delta[1], ind_delta[2])
        llik.prop = logPostDelta(delta_prop)
        llik.cur = logPostDelta(delta_cur)
        logr1 = min(0, (llik.prop-llik.cur+
                          dbeta(delta_cur, shape1 = ind_delta[1], shape2 = ind_delta[2], log = TRUE) -
                          dbeta(delta_prop, shape1 = ind_delta[1], shape2 = ind_delta[2], log = TRUE)))
      }
      if(prop_delta == "RW"){
        lgdelta_cur = log(delta_cur/(1-delta_cur))
        lgdelta_prop = rnorm(1, mean = lgdelta_cur, sd = rw_delta)
        delta_prop = exp(lgdelta_prop)/(1+exp(lgdelta_prop))
        llik.prop = logPostDelta(delta_prop)
        llik.cur = logPostDelta(delta_cur)
        logr1 = min(0, (llik.prop-llik.cur+
                          log(delta_prop)+log(1-delta_prop)-log(delta_cur)-log(1-delta_cur)))
      }
      if(runif(1) <= exp(logr1)){
        delta_cur = delta_prop; counter1 = counter1+1
      }
      if( i > burnin & (i-burnin) %% thin==0) {
        delta_chain[(i-burnin)/thin] = delta_cur
      }
      # Update beta
      if(prop_beta == "NORMAL"){
        beta_prop = rnorm(k, mean = c(t(beta_cur)), sd = rw_n_beta)
        logr2 = min(0, logPostBeta(beta_prop)-logPostBeta(beta_cur))
      }
      if(prop_beta == "UNIF"){
        beta_prop = runif(k, min = beta_cur-rw_u_beta, max = beta_cur+rw_u_beta)
        logr2 = min(0, logPostBeta(beta_prop)-logPostBeta(beta_cur))
      }
      if(runif(1) <= exp(logr2)){
        beta_cur = beta_prop; counter2 = counter2+1
      }
      if( i > burnin & (i-burnin) %% thin==0) {
        beta_chain[(i-burnin)/thin, ] = beta_cur
      }
    }
    return(list(sample = cbind(beta_chain, delta_chain), 
                accept = c(counter1/n_mcmc, counter2/n_mcmc)))
  }
  
  outcome_NPP <- lr_NPP_mcmc(logC_seq=list(delta=delta, logC=logC))
  
  accept <- outcome_NPP$accept
  meanout <- colMeans(outcome_NPP$sample[,1:(k+1)])
  sdout <- apply(outcome_NPP$sample[,1:(k+1)],2,sd)
  mdout <- apply(outcome_NPP$sample[,1:(k+1)],2,median)
  mod <- KernSmooth::bkde(outcome_NPP$sample[,(k+1)])
  modedelta <- mod$x[which(mod$y == max(mod$y))]
  return(c(accept,meanout,sdout,mdout,modedelta))
}
