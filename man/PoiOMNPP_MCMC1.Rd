\name{PoiOMNPP_MCMC1}
\alias{PoiOMNPP_MCMC1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC Sampling for Poisson Population of multiple ordered historical data using Normalized Power Prior}

\description{
Multiple ordered historical data are incorporated together.
Conduct posterior sampling for Poisson population with normalized power prior.
For the power parameter \eqn{\gamma}, a Metropolis-Hastings algorithm with independence proposal is used.
For the model parameter \eqn{\lambda}, Gibbs sampling is used.
}

\usage{
PoiOMNPP_MCMC1(n0,n,prior_gamma,prior_lambda, gamma_ind_prop,
               gamma_ini,nsample,burnin,thin)
}

\arguments{
\item{n0}{a natural number vector : number of successes in historical data.}
\item{n}{a natural number : number of successes in the current data.}
\item{prior_gamma}{a vector of the hyperparameters in the prior distribution \eqn{Dirichlet(\alpha_1, \alpha_2, ... ,\alpha_K)} for \eqn{\gamma}.}
\item{prior_lambda}{a vector of the hyperparameters in the prior distribution \eqn{Gamma(\alpha, \beta)} for \eqn{\lambda}.}
\item{gamma_ind_prop}{a vector of the hyperparameters in the proposal distribution \eqn{Dirichlet(\alpha_1, \alpha_2, ... ,\alpha_K)} for \eqn{\gamma}.}
\item{gamma_ini}{the initial value of \eqn{\gamma} in MCMC sampling.}
\item{nsample}{specifies the number of posterior samples in the output.}
\item{burnin}{the number of burn-ins. The output will only show MCMC samples after bunrin.}
\item{thin}{the thinning parameter in MCMC sampling.}
}

\value{
A list of class "NPP" with three elements:
\item{acceptrate}{the acceptance rate in MCMC sampling for \eqn{\gamma} using Metropolis-Hastings algorithm.}
\item{lambda}{posterior of the model parameter \eqn{\lambda}.}
\item{delta}{posterior of the power parameter \eqn{\delta}. It is equal to the cumulative sum of \eqn{\gamma}}
}

\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\gamma}.
The normalized power prior distribution is
\deqn{\frac{\pi_0(\gamma)\pi_0(\lambda)\prod_{k=1}^{K}L(\lambda|D_{0k})^{(\sum_{i=1}^{k}\gamma_i)}}{\int \pi_0(\lambda)\prod_{k=1}^{K}L(\lambda|D_{0k})^{(\sum_{i=1}^{k}\gamma_i)}d\lambda }.}

Here \eqn{\pi_0(\gamma)} and \eqn{\pi_0(\lambda)} are the initial prior distributions of \eqn{\gamma} and \eqn{\lambda}, respectively. \eqn{L(\lambda|D_{0k})} is the likelihood function of historical data \eqn{D_{0k}}, and \eqn{\sum_{i=1}^{k}\gamma_i} is the corresponding power parameter.
}

\examples{
PoiOMNPP_MCMC1(n0=c(0,3,5),n=3,prior_gamma=c(1/2,1/2,1/2,1/2), prior_lambda=c(1,1/10),
               gamma_ind_prop=rep(1,4),gamma_ini=NULL, nsample = 2000, burnin = 500, thin = 2)
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
         \code{\link{PoiMNPP_MCMC2}};
         \code{\link{PoiOMNPP_MCMC2}}
}

\keyword{multiple historical data normalized power prior}
