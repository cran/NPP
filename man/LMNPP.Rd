\name{LMNPP_MCMC}
\alias{LMNPP_MCMC}
\title{MCMC Sampling for Normal Linear Model using Normalized Power Prior}


\description{
Conduct posterior sampling for normal linear model with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either
independence proposal, or a random walk proposal on its logit scale is used.
For the regression parameter \eqn{\beta} and \eqn{\sigma^2}, Gibbs sampling is used.
}


\usage{
LMNPP_MCMC(y.Cur, y.Hist, x.Cur = NULL, x.Hist = NULL,
           prior = list(a = 1.5, b = 0, mu0 = 0,
                   Rinv = matrix(1, nrow = 1), delta.alpha = 1, delta.beta = 1),
           MCMCmethod = 'IND', rw.logit.delta = 0.1,
           ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000,
           control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
}

\arguments{
\item{y.Cur}{a vector of individual level of the response y in current data.}

\item{y.Hist}{a vector of individual level of the response y in historical data.}


\item{x.Cur}{a vector or matrix or data frame of covariate observed in the current data.
             If more than 1 covariate available, the number of rows is equal to the number of
             observations.}

\item{x.Hist}{a vector or matrix or data frame of covariate observed in the historical data.
             If more than 1 covariate available, the number of rows is equal to the number of
             observations.}


\item{prior}{a list of the hyperparameters in the prior for model parameters
             \eqn{(\beta, \sigma^2)} and \eqn{\delta}.
             The form of the prior for model parameter \eqn{(\beta, \sigma^2)} is in the section "Details".

     \code{a} a positive hyperparameter for prior on model parameters. It is the power \eqn{a} in formula \eqn{(1/\sigma^2)^a};
     See details.

     \code{b} equals 0 if a flat prior is used for \eqn{\beta}. Equals 1 if a normal prior is used for \eqn{\beta}; See details.

     \code{mu0} a vector of the mean for prior \eqn{\beta|\sigma^2}. Only applicable if \code{b = 1}.

     \code{Rinv} inverse of the matrix \eqn{R}. The covariance matrix of the prior for \eqn{\beta|\sigma^2} is \eqn{\sigma^2 R^{-1}}.

     \code{delta.alpha} is the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

     \code{delta.beta} is the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

}


\item{MCMCmethod}{sampling method for \eqn{\delta} in MCMC. It can be either 'IND' for independence proposal; or 'RW' for
random walk proposal on logit scale. }

\item{rw.logit.delta}{the stepsize(variance of the normal distribution) for the random walk proposal of logit \eqn{\delta}.
                      Only applicable if \code{MCMCmethod = 'RW'. } }

\item{ind.delta.alpha}{specifies the first parameter \eqn{\alpha} when
independent proposal \eqn{Beta(\alpha, \beta)} for \eqn{\delta} is used. Only applicable if \code{MCMCmethod = 'IND'}}


\item{ind.delta.beta}{specifies the first parameter \eqn{\beta} when
independent proposal \eqn{Beta(\alpha, \beta)} for \eqn{\delta} is used. Only applicable if \code{MCMCmethod = 'IND'}}


\item{nsample}{specifies the number of posterior samples in the output.}

\item{control.mcmc}{a list of three elements used in posterior sampling.

       \code{delta.ini} is the initial value of \eqn{\delta} in MCMC sampling.

       \code{burnin} is the number of burn-ins. The output will only show MCMC samples after bunrin.

       \code{thin} is the thinning parameter in MCMC sampling.
       }
}


\value{
A list of class "NPP" with five elements:
\item{beta}{posterior of the model parameter \eqn{\beta} in vector or matrix form.}

\item{sigmasq}{posterior of the model parameter \eqn{\sigma^2}.}

\item{delta}{posterior of the power parameter \eqn{\delta}.}

\item{acceptance}{the acceptance rate in MCMC sampling for \eqn{\delta} using Metropolis-Hastings algorithm.}

\item{DIC}{the deviance information criteria for model diagnostics.}

}




\details{
If \eqn{b = 1}, prior for \eqn{(\beta, \sigma)} is \eqn{(1/\sigma^2)^a * N(mu0, \sigma^2 R^{-1})}, which includes the g-prior.
If \eqn{b = 0}, prior for \eqn{(\beta, \sigma)} is \eqn{(1/\sigma^2)^a}.
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate when sampling \eqn{\delta}, and
the deviance information criteria.
}



\examples{
set.seed(123)
x1 = runif(100, min = 0, max = 10)
x0 = runif(100, min = 0, max = 1)
y1 = 10+ 2*x1 + rnorm(100, mean = 0, sd = 1)
y0 = 10+ 1.5*x0 + rnorm(100, mean = 0, sd = 1)

RegPost = LMNPP_MCMC(y.Cur = y1, y.Hist = y0, x.Cur = x1, x.Hist = x0,
                     prior = list(a = 1.5, b = 0, mu0 = c(0, 0),
                                  Rinv = diag(100, nrow = 2),
                     delta.alpha = 1, delta.beta = 1), MCMCmethod = 'IND',
                     ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000,
                     control.mcmc = list(delta.ini = NULL,
                                         burnin = 2000, thin = 2))
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}


\references{

Ibrahim, J.G., Chen, M.-H., Gwon, Y. and Chen, F. (2015).
The Power Prior: Theory and Applications.
\emph{Statistics in Medicine  34:3724-3749.}


Duan, Y., Ye, K. and Smith, E.P. (2006).
Evaluating Water Quality: Using Power Priors to Incorporate Historical Information.
\emph{Environmetrics  17:95-106.}


Berger, J.O. and Bernardo, J.M. (1992).
On the development of reference priors.
\emph{Bayesian Statistics 4: Proceedings of
the Fourth Valencia International Meeting, Bernardo, J.M, Berger,
J.O., Dawid, A.P. and Smith, A.F.M. eds., 35-60, Clarendon Press:Oxford.}


Jeffreys, H. (1946).
An Invariant Form for the Prior Probability in Estimation Problems.
\emph{Proceedings of the Royal Statistical Society of London, Series A  186:453-461.}

}


\seealso{
         \code{\link{BerNPP_MCMC}};
         \code{\link{MultinomialNPP_MCMC}};
         \code{\link{PoissonNPP_MCMC}};
         \code{\link{NormalNPP_MCMC}}
}

\keyword{normalized power prior}
