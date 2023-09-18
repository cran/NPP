\name{PoiMNPP_MCMC1}
\alias{PoiMNPP_MCMC1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC Sampling for Poisson Population using Normalized Power Prior with Multiple Historical Data}

\description{
This function incorporates multiple sets of historical data for posterior sampling in a Poisson population using a normalized power prior. The power parameter \eqn{\delta} uses a Metropolis-Hastings algorithm, which can be either an independence proposal or a random walk proposal on its logit scale. For the model parameter \eqn{\lambda}, Gibbs sampling is employed.
}

\usage{
PoiMNPP_MCMC1(n0, n, prior_lambda, prop_delta, prior_delta_alpha,
              prior_delta_beta, rw_delta, delta_ini, nsample, burnin, thin)
}

\arguments{
\item{n0}{A vector of natural numbers: number of successes in historical data.}

\item{n}{A natural number: number of successes in the current data.}

\item{prior_lambda}{A vector of hyperparameters for the prior distribution \eqn{Gamma(\alpha, \beta)} of \eqn{\lambda}.}

\item{prop_delta}{The class of proposal distribution for \eqn{\delta}.}

\item{prior_delta_alpha}{A vector of hyperparameter \eqn{\alpha} for the prior distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}

\item{prior_delta_beta}{A vector of hyperparameter \eqn{\beta} for the prior distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}

\item{rw_delta}{The stepsize (variance of the normal distribution) for the random walk proposal of logit \eqn{\delta}. This is only applicable if \code{prop_delta = 'RW'}.}

\item{delta_ini}{The initial value for \eqn{\delta} in MCMC sampling.}

\item{nsample}{Specifies the number of posterior samples in the output.}

\item{burnin}{The number of burn-ins. Only the MCMC samples after this burn-in will be shown in the output.}

\item{thin}{The thinning parameter used in MCMC sampling.}
}

\value{
A list of class "NPP" comprising:
\item{acceptrate}{The acceptance rate in MCMC sampling for \eqn{\delta} using the Metropolis-Hastings algorithm.}
\item{lambda}{Posterior samples of the model parameter \eqn{\lambda}.}
\item{delta}{Posterior samples of the power parameter \eqn{\delta}.}
}

\details{
The function returns posteriors for both the model and power parameters, as well as the acceptance rate for sampling \eqn{\delta}. The normalized power prior distribution is given by:
\deqn{\frac{\pi_0(\delta)\pi_0(\lambda)\prod_{k=1}^{K}L(\lambda|D_{0k})^{\delta_{k}}}{\int \pi_0(\lambda)\prod_{k=1}^{K}L(\lambda|D_{0k})^{\delta_{k}} d\lambda}.}

Here, \eqn{\pi_0(\delta)} and \eqn{\pi_0(\lambda)} are the initial prior distributions for \eqn{\delta} and \eqn{\lambda}, respectively. \eqn{L(\lambda|D_{0k})} is the likelihood function based on historical data \eqn{D_{0k}}, with \eqn{\delta_k} being its corresponding power parameter.
}

\examples{
PoiMNPP_MCMC1(n0 = c(0, 3, 5), n = 3, prior_lambda = c(1, 1/10), prop_delta = "IND",
              prior_delta_alpha = c(1, 1, 1), prior_delta_beta = c(1, 1, 1),
              rw_delta = 0.1, delta_ini = NULL, nsample = 2000, burnin = 500, thin = 2)
}

\author{
Qiang Zhang \email{zqzjf0408@163.com}
}

\references{
Ibrahim, J.G., Chen, M.-H., Gwon, Y., and Chen, F. (2015).
The Power Prior: Theory and Applications.
\emph{Statistics in Medicine  34:3724-3749.}

Duan, Y., Ye, K., and Smith, E.P. (2006).
Evaluating Water Quality: Using Power Priors to Incorporate Historical Information.
\emph{Environmetrics  17:95-106.}
}

\seealso{
\code{\link{PoiMNPP_MCMC2}},
\code{\link{PoiOMNPP_MCMC1}},
\code{\link{PoiOMNPP_MCMC2}}
}

\keyword{multiple historical data}
\keyword{normalized power prior}
