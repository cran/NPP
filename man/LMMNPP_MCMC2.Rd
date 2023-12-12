\name{LMMNPP_MCMC2}
\alias{LMMNPP_MCMC2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC Sampling for Linear Regression Model of multiple historical data using
Normalized Power Prior
}

\description{
Multiple historical data are combined individually.
The NPP of multiple historical data is the product of the NPP of each historical data.
Conduct posterior sampling for Linear Regression Model with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either
independence proposal, or a random walk proposal on its logit scale is used.
For the model parameters \eqn{(\beta, \sigma^2)}, Gibbs sampling is used.
}

\usage{
LMMNPP_MCMC2(D0, X, Y, a0, b, mu0, R, delta_ini, prop_delta,
             prior_delta_alpha, prior_delta_beta, prop_delta_alpha,
             prop_delta_beta, rw_delta, nsample, burnin, thin)
}

\arguments{
\item{D0}{a list of \eqn{k} elements representing \eqn{k} historical data, where the
          \eqn{i^{th}} element corresponds to the \eqn{i^{th}} historical data named
          as ``D0i''.}
\item{X}{a vector or matrix or data frame of covariate observed in the current data.
         If more than 1 covariate available, the number of rows is equal to the
         number of observations.}
\item{Y}{a vector of individual level of the response y in the current data.}
\item{a0}{a positive shape parameter for inverse-gamma prior on model parameter                          \eqn{\sigma^2}.}
\item{b}{a positive scale parameter for inverse-gamma prior on model parameter \eqn{\sigma^2}.}
\item{mu0}{a vector of the mean for prior \eqn{\beta|\sigma^2}.}
\item{R}{a inverse matrix of the covariance matrix for prior \eqn{\beta|\sigma^2}.}
\item{delta_ini}{the initial value of \eqn{\delta} in MCMC sampling.}
\item{prop_delta}{the class of proposal distribution for \eqn{\delta}.}
\item{prior_delta_alpha}{a vector of the hyperparameter \eqn{\alpha} in the prior distribution                          \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}
\item{prior_delta_beta}{a vector of the hyperparameter \eqn{\beta} in the prior distribution                           \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}
\item{prop_delta_alpha}{a vector of the hyperparameter \eqn{\alpha} in the proposal                                    distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}
\item{prop_delta_beta}{a vector of the hyperparameter \eqn{\beta} in the proposal distribution                        \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}
\item{rw_delta}{the stepsize(variance of the normal distribution) for the random walk proposal
                of logit \eqn{\delta}. Only applicable if prop_delta = 'RW'.}
\item{nsample}{specifies the number of posterior samples in the output.}
\item{burnin}{the number of burn-ins. The output will only show MCMC samples after bunrin.}
\item{thin}{the thinning parameter in MCMC sampling.}
}

\details{
The outputs include posteriors of the model parameters and power parameter,
acceptance rate in sampling \eqn{\delta}.
Let \eqn{\theta}=\eqn{(\beta, \sigma^2)}, the normalized power prior distribution is
\deqn{\pi_0(\delta)\prod_{k=1}^{K}\frac{\pi_0(\theta)L(\theta|D_{0k})^{\delta_k}}{\int \pi_0(\theta)L(\theta|D_{0k})^{\delta_k} \,d\theta}.}

Here \eqn{\pi_0(\delta)} and \eqn{\pi_0(\theta)} are the initial prior distributions of \eqn{\delta} and \eqn{\theta}, respectively. \eqn{L(\theta|D_{0k})} is the likelihood function of historical data \eqn{D_{0k}}, and \eqn{\delta_k} is the corresponding power parameter.
}

\value{
A list of class "NPP" with four elements:
\item{acceptrate}{the acceptance rate in MCMC sampling for \eqn{\delta} using
                  Metropolis-Hastings algorithm.}
\item{beta}{posterior of the model parameter \eqn{\beta} in vector or matrix form.}
\item{sigma}{posterior of the model parameter \eqn{\sigma^2}.}
\item{delta}{posterior of the power parameter \eqn{\delta}.}
}

\examples{
\dontrun{
set.seed(1234)
sigsq0 = 1

n01 = 100
theta01 = c(0, 1, 1)
X01 = cbind(1, rnorm(n01, mean=0, sd=1), runif(n01, min=-1, max=1))
Y01 = X01\%*\%as.vector(theta01) + rnorm(n01, mean=0, sd=sqrt(sigsq0))
D01 = cbind(X01, Y01)

n02 = 70
theta02 = c(0, 2, 3)
X02 = cbind(1, rnorm(n02, mean=0, sd=1), runif(n02, min=-1, max=1))
Y02 = X02\%*\%as.vector(theta02) + rnorm(n02, mean=0, sd=sqrt(sigsq0))
D02 = cbind(X02, Y02)

n03 = 50
theta03 = c(0, 3, 5)
X03 = cbind(1, rnorm(n03, mean=0, sd=1), runif(n03, min=-1, max=1))
Y03 = X03\%*\%as.vector(theta03) + rnorm(n03, mean=0, sd=sqrt(sigsq0))
D03 = cbind(X03, Y03)

D0 = list(D01, D02, D03)
n0 = c(n01, n02, n03)

n = 100
theta = c(0, 3, 5)
X = cbind(1, rnorm(n, mean=0, sd=1), runif(n, min=-1, max=1))
Y = X\%*\%as.vector(theta) + rnorm(n, mean=0, sd=sqrt(sigsq0))

LMMNPP_MCMC2(D0=D0, X=X, Y=Y, a0=2, b=2, mu0=c(0,0,0), R=diag(c(1/64,1/64,1/64)),
             delta_ini=NULL, prior_delta_alpha=c(1,1,1), prior_delta_beta=c(1,1,1),
             prop_delta_alpha=c(1,1,1), prop_delta_beta=c(1,1,1),
             prop_delta="RW", rw_delta=0.9, nsample=5000, burnin=1000, thin=5)
}
}

\author{
Qiang Zhang \email{zqzjf0408@163.com}
}

\references{
Ibrahim, J.G., Chen, M.-H., Gwon, Y. and Chen, F. (2015).
The Power Prior: Theory and Applications.
\emph{Statistics in Medicine  34:3724-3749.}

Duan, Y., Ye, K. and Smith, E.P. (2006).
Evaluating Water Quality: Using Power Priors to Incorporate Historical Information.
\emph{Environmetrics  17:95-106.}
}

\seealso{
         \code{\link{LMMNPP_MCMC1}};
         \code{\link{LMOMNPP_MCMC1}};
         \code{\link{LMOMNPP_MCMC2}}
}

\keyword{multiple historical data normalized power prior}
