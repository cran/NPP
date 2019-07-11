\name{PoissonNPP_MCMC}
\alias{PoissonNPP_MCMC}
\title{MCMC Sampling for Bernoulli Population using Normalized Power Prior}


\description{
Conduct posterior sampling for Poisson population with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either
independence proposal, or a random walk proposal on its logit scale is used.
For the model parameter \eqn{\lambda}, Gibbs sampling is used.
}


\usage{
PoissonNPP_MCMC(Data.Cur, Data.Hist,
                CompStat = list(n0 = NULL, mean0 = NULL, n1 = NULL, mean1 = NULL),
                prior = list(lambda.shape = 1/2, lambda.scale = 100,
                             delta.alpha = 1, delta.beta = 1),
                MCMCmethod = 'IND', rw.logit.delta = 0.1,
                ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000,
                control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
}

\arguments{
\item{Data.Cur}{a non-negative integer vector of each observed current data.}

\item{Data.Hist}{a non-negative integer vector of each observed historical data.}

\item{CompStat}{a list of four elements that represents the
                "compatibility(sufficient) statistics" for \eqn{\lambda}.
                Default is \code{NULL} so the fitting will be based on
                the data. If the \code{CompStat} is provided then the inputs in
                \code{Data.Cur} and \code{Data.Hist} will be ignored.

     \code{n0} is the number of observations in the historical data.

     \code{mean0} is the sample mean of the historical data.

     \code{n1} is the number of observations in the current data.

     \code{mean1} is the sample mean of the current data.
}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{\lambda} and \eqn{\delta}.
             A Gamma distribution is used as the prior of \eqn{\lambda}, and a Beta distribution is used as the prior of \eqn{\delta}.

     \code{lambda.shape} is the shape (hyper)parameter in the prior distribution \eqn{Gamma(shape, scale)} for \eqn{\lambda}.

     \code{lambda.scale} is the scale (hyper)parameter in the prior distribution \eqn{Gamma(shape, scale)} for \eqn{\lambda}.

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
\item{lambda}{posterior of the model parameter \eqn{\lambda}.}

\item{delta}{posterior of the power parameter \eqn{\delta}.}

\item{acceptance}{the acceptance rate in MCMC sampling for \eqn{\delta} using Metropolis-Hastings algorithm.}

\item{DIC}{the deviance information criteria for model diagnostics.}

}





\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\delta}, and
the deviance information criteria.
}



\examples{
set.seed(1234)
DataHist <- rpois(n = 100, lambda = 49)
set.seed(12345)
DataCur <- rpois(n = 100, lambda = 49)

PoissonNPP_MCMC(Data.Cur = DataCur, Data.Hist = DataHist,
                CompStat = list(n0 = 20, mean0 = 10, n1 = 30, mean1 = 11),
                prior = list(lambda.shape = 1/2, lambda.scale = 100,
                             delta.alpha = 1, delta.beta = 1),
                MCMCmethod = 'RW', rw.logit.delta = 1,
                ind.delta.alpha= 1, ind.delta.beta= 1,nsample = 10000,
                control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 1))
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
         \code{\link{BerNPP_MCMC}};
}

\keyword{normalized power prior}
