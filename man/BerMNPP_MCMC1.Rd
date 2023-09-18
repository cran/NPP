\name{BerMNPP_MCMC1}
\alias{BerMNPP_MCMC1}
\title{MCMC Sampling for Bernoulli Population with Multiple Historical Data using Normalized Power Prior}

\description{
  Incorporate multiple historical data sets for posterior sampling of a Bernoulli population using the normalized power prior. The Metropolis-Hastings algorithm, with either an independence proposal or a random walk proposal on the logit scale, is applied for the power parameter \eqn{\delta}. Gibbs sampling is utilized for the model parameter \eqn{p}.
}

\usage{
  BerMNPP_MCMC1(n0, y0, n, y, prior_p, prior_delta_alpha,
                prior_delta_beta, prop_delta_alpha, prop_delta_beta,
                delta_ini, prop_delta, rw_delta, nsample, burnin, thin)
}

\arguments{
\item{n0}{A non-negative integer vector representing the number of trials in historical data.}

\item{y0}{A non-negative integer vector denoting the number of successes in historical data.}

\item{n}{A non-negative integer indicating the number of trials in the current data.}

\item{y}{A non-negative integer for the number of successes in the current data.}

\item{prior_p}{a vector of the hyperparameters in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{p}.}

\item{prior_delta_alpha}{a vector of the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}

\item{prior_delta_beta}{a vector of the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}

\item{prop_delta_alpha}{a vector of the hyperparameter \eqn{\alpha} in the proposal distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}

\item{prop_delta_beta}{a vector of the hyperparameter \eqn{\beta} in the proposal distribution \eqn{Beta(\alpha, \beta)} for each \eqn{\delta}.}

\item{delta_ini}{the initial value of \eqn{\delta} in MCMC sampling.}

\item{prop_delta}{the class of proposal distribution for \eqn{\delta}.}

\item{rw_delta}{the stepsize(variance of the normal distribution) for the random walk proposal of logit \eqn{\delta}. Only applicable if prop_delta = 'RW'.}

\item{nsample}{specifies the number of posterior samples in the output.}

\item{burnin}{the number of burn-ins. The output will only show MCMC samples after bunrin.}

\item{thin}{the thinning parameter in MCMC sampling.}
}

\value{
  A list of class "NPP" comprising:
    \item{acceptrate}{Acceptance rate in MCMC sampling for \eqn{\delta} via the Metropolis-Hastings algorithm.}
  \item{p}{Posterior distribution of the model parameter \eqn{p}.}
  \item{delta}{Posterior distribution of the power parameter \eqn{\delta}.}
}

\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\delta}.
The normalized power prior distribution is
\deqn{\frac{\pi_0(\delta)\pi_0(\theta)\prod_{k=1}^{K}L(\theta|D_{0k})^{\delta_{k}}}{\int \pi_0(\theta)\prod_{k=1}^{K}L(\theta|D_{0k})^{\delta_{k}} d\theta}.}

Here \eqn{\pi_0(\delta)} and \eqn{\pi_0(\theta)} are the initial prior distributions of \eqn{\delta} and \eqn{\theta}, respectively. \eqn{L(\theta|D_{0k})} is the likelihood function of historical data \eqn{D_{0k}}, and \eqn{\delta_k} is the corresponding power parameter.
}

\examples{
BerMNPP_MCMC1(n0 = c(275, 287), y0 = c(92, 125), n = 39, y = 17,
              prior_p = c(1/2,1/2), prior_delta_alpha = c(1/2,1/2),
              prior_delta_beta = c(1/2,1/2),
              prop_delta_alpha = c(1,1)/2, prop_delta_beta = c(1,1)/2,
              delta_ini = NULL, prop_delta = "IND",
              nsample = 2000, burnin = 500, thin = 2)
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
         \code{\link{BerMNPP_MCMC2}};
         \code{\link{BerOMNPP_MCMC1}};
         \code{\link{BerOMNPP_MCMC2}}
}

\keyword{multiple historical data normalized power prior}
