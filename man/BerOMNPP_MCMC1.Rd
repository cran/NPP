\name{BerOMNPP_MCMC1}
\alias{BerOMNPP_MCMC1}

\title{MCMC Sampling for Bernoulli Population of multiple ordered historical data using Normalized Power Prior}

\description{
Multiple ordered historical data are incorporated together.
Conduct posterior sampling for Bernoulli population with normalized power prior.
For the power parameter \eqn{\gamma}, a Metropolis-Hastings algorithm with independence proposal is used.
For the model parameter \eqn{p}, Gibbs sampling is used.
}

\usage{
    BerOMNPP_MCMC1(n0, y0, n, y, prior_gamma, prior_p, gamma_ind_prop,
                   gamma_ini, nsample, burnin, thin, adjust = FALSE)
}

\arguments{
\item{n0}{a non-negative integer vector: number of trials in historical data.}

\item{y0}{a non-negative integer vector: number of successes in historical data.}

\item{n}{a non-negative integer: number of trials in the current data.}

\item{y}{a non-negative integer: number of successes in the current data.}

\item{prior_gamma}{a vector of the hyperparameters in the prior distribution \eqn{Dirichlet(\alpha_1, \alpha_2, ... ,\alpha_K)} for \eqn{\gamma}.}

\item{prior_p}{a vector of the hyperparameters in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{p}.}

\item{gamma_ind_prop}{a vector of the hyperparameters in the proposal distribution \eqn{Dirichlet(\alpha_1, \alpha_2, ... ,\alpha_K)} for \eqn{\gamma}.}

\item{gamma_ini}{the initial value of \eqn{\gamma} in MCMC sampling.}

\item{nsample}{specifies the number of posterior samples in the output.}

\item{burnin}{the number of burn-ins. The output will only show MCMC samples after burnin.}

\item{thin}{the thinning parameter in MCMC sampling.}

\item{adjust}{Logical, indicating whether or not to adjust the parameters of the proposal distribution.}
}

\value{
A list of class "NPP" with three elements:
\item{acceptrate}{the acceptance rate in MCMC sampling for \eqn{\gamma} using Metropolis-Hastings algorithm.}
\item{p}{posterior of the model parameter \eqn{p}.}
\item{delta}{posterior of the power parameter \eqn{\delta}. It is equal to the cumulative sum of \eqn{\gamma}.}
}

\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\gamma}.
The normalized power prior distribution is given by:
\deqn{\frac{\pi_0(\gamma)\pi_0(\theta)\prod_{k=1}^{K}L(\theta|D_{0k})^{(\sum_{i=1}^{k}\gamma_i)}}{\int \pi_0(\theta)\prod_{k=1}^{K}L(\theta|D_{0k})^{(\sum_{i=1}^{k}\gamma_i)}d\theta}.}
Here, \eqn{\pi_0(\gamma)} and \eqn{\pi_0(\theta)} are the initial prior distributions of \eqn{\gamma} and \eqn{\theta}, respectively. \eqn{L(\theta|D_{0k})} is the likelihood function of historical data \eqn{D_{0k}}, and \eqn{\sum_{i=1}^{k}\gamma_i} is the corresponding power parameter.
}

\examples{
BerOMNPP_MCMC1(n0 = c(275, 287), y0 = c(92, 125), n = 39, y = 17, prior_gamma=c(1,1,1)/3,
               prior_p=c(1/2,1/2), gamma_ind_prop=rep(1,3)/2, gamma_ini=NULL,
               nsample = 2000, burnin = 500, thin = 2, adjust = FALSE)
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
\code{\link{BerMNPP_MCMC1}},
\code{\link{BerMNPP_MCMC2}},
\code{\link{BerOMNPP_MCMC2}}
}

\keyword{multiple historical data normalized power prior}
