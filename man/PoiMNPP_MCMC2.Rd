\name{PoiMNPP_MCMC2}
\alias{PoiMNPP_MCMC2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC Sampling for Poisson Population of multiple historical data using Normalized Power Prior}

\description{
Multiple historical data are combined individually.
Conduct posterior sampling for Poisson population with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either independence proposal, or a random walk proposal on its logit scale is used.
For the model parameter \eqn{\lambda}, Gibbs sampling is used.
}

\usage{
PoiMNPP_MCMC2(n0,n,prior_lambda,prop_delta,prior_delta_alpha,
                  prior_delta_beta,rw_delta, delta_ini,nsample,burnin,thin)
}

\arguments{
\item{n0}{a natural number vector : number of successes in historical data.}
\item{n}{a natural number : number of successes in the current data.}
\item{prior_lambda}{a vector of the hyperparameters in the prior distribution \eqn{Gamma(\alpha, \beta)} for \eqn{\lambda}.}
\item{prop_delta}{the class of proposal distribution for \eqn{\delta}.}
\item{prior_delta_alpha}{a vector of the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}
\item{prior_delta_beta}{a vector of the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}
\item{rw_delta}{the stepsize(variance of the normal distribution) for the random walk proposal of logit \eqn{\delta}. Only applicable if prop_delta = 'RW'.}
\item{delta_ini}{the initial value of \eqn{\delta} in MCMC sampling.}
\item{nsample}{specifies the number of posterior samples in the output.}
\item{burnin}{the number of burn-ins. The output will only show MCMC samples after bunrin.}
\item{thin}{the thinning parameter in MCMC sampling.}
}

\value{
A list of class "NPP" with three elements:
\item{acceptrate}{the acceptance rate in MCMC sampling for \eqn{\delta} using Metropolis-Hastings algorithm.}
\item{lambda}{posterior of the model parameter \eqn{\lambda}.}
\item{delta}{posterior of the power parameter \eqn{\delta}.}
}

\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\delta}.
The normalized power prior distribution is
\deqn{\pi_0(\delta)\prod_{k=1}^{K}\frac{\pi_0(\lambda)L(\lambda|D_{0k})^{\delta_{k}}}{\int \pi_0(\lambda)L(\lambda|D_{0k})^{\delta_{k}} d\lambda}.}

Here \eqn{\pi_0(\delta)} and \eqn{\pi_0(\lambda)} are the initial prior distributions of \eqn{\delta} and \eqn{\lambda}, respectively. \eqn{L(\lambda|D_{0k})} is the likelihood function of historical data \eqn{D_{0k}}, and \eqn{\delta_k} is the corresponding power parameter.
}

\examples{
PoiMNPP_MCMC2(n0=c(0,3,5),n=3,prior_lambda=c(1,1/10),prop_delta="IND",
              prior_delta_alpha=c(1,1,1), prior_delta_beta=c(1,1,1), rw_delta=0.1,
              delta_ini=NULL, nsample = 2000, burnin = 500, thin = 2)
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
         \code{\link{PoiMNPP_MCMC1}};
         \code{\link{PoiOMNPP_MCMC1}};
         \code{\link{PoiOMNPP_MCMC2}}
}

\keyword{multiple historical data normalized power prior}
