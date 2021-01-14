\name{BerNPP_MCMC}
\alias{BerNPP_MCMC}
\title{MCMC Sampling for Bernoulli Population using Normalized Power Prior}


\description{
Conduct posterior sampling for Bernoulli population with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either
independence proposal, or a random walk proposal on its logit scale is used.
For the model parameter \eqn{p}, Gibbs sampling is used.
}


\usage{
BerNPP_MCMC(Data.Cur = c(100, 50), Data.Hist = c(100, 50),
            CompStat = list(n0 = NULL, y0 = NULL, n1 = NULL, y1 = NULL),
            prior = list(p.alpha = 1, p.beta = 1, delta.alpha = 1, delta.beta = 1),
            MCMCmethod = 'IND', rw.logit.delta = 0.1,
            ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000,
            control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
}

\arguments{
\item{Data.Cur}{a non-negative integer vector of two elements: c(number of trials, number of successes) in the current data.}

\item{Data.Hist}{a non-negative integer vector of two elements: c(number of trials, number of successes) in the historical data.}

\item{CompStat}{a list of four elements that represents the
                "compatibility(sufficient) statistics" for \eqn{p}.
                Default is \code{NULL} so the fitting will be based on
                the data. If the \code{CompStat} is provided then the inputs in
                \code{Data.Cur} and \code{Data.Hist} will be ignored.
                Note: in Bernoulli population providing \code{CompStat} is equivalent to provide the data summary as in
                \code{Data.Cur} and \code{Data.Cur}.

     \code{n0} is the number of trials in the historical data.

     \code{y0} is the number of successes in the historical data.

     \code{n1} is the number of trials in the current data.

     \code{y1} is the number of successes in the current data.
}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{p} and \eqn{\delta}.

     \code{p.alpha} is the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{p}.

     \code{p.beta} is the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{p}.

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
A list of class "NPP" with four elements:
\item{p}{posterior of the model parameter \eqn{p}.}

\item{delta}{posterior of the power parameter \eqn{\delta}.}

\item{acceptance}{the acceptance rate in MCMC sampling for \eqn{\delta} using Metropolis-Hastings algorithm. }

\item{DIC}{the deviance information criteria for model diagnostics.}

}





\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\delta}, and
the deviance information criteria.
}



\examples{
BerNPP_MCMC(Data.Cur = c(493, 473), Data.Hist = c(680, 669),
            prior = list(p.alpha = 0.5, p.beta = 0.5, delta.alpha = 1, delta.beta = 1),
            MCMCmethod = 'RW', rw.logit.delta = 1, nsample = 5000,
            control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
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

}


\seealso{
         \code{\link{MultinomialNPP_MCMC}};
         \code{\link{NormalNPP_MCMC}};
         \code{\link{PoissonNPP_MCMC}}
}

\keyword{normalized power prior}
