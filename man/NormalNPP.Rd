\name{NormalNPP_MCMC}
\alias{NormalNPP_MCMC}
\title{MCMC Sampling for Normal Population using Normalized Power Prior}


\description{
Conduct posterior sampling for normal population with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either
independence proposal, or a random walk proposal on its logit scale is used.
For the model parameter \eqn{\mu} and \eqn{\sigma^2}, Gibbs sampling is used.
}


\usage{
NormalNPP_MCMC(Data.Cur, Data.Hist,
               CompStat = list(n0 = NULL, mean0 = NULL, var0 = NULL,
                               n1 = NULL, mean1 = NULL, var1 = NULL),
               prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1),
               MCMCmethod = 'IND', rw.logit.delta = 0.1,
               ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 5000,
               control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
}

\arguments{
\item{Data.Cur}{a vector of individual level current data.}

\item{Data.Hist}{a vector of individual level historical data.}

\item{CompStat}{a list of six elements(scalar) that represents the
                "compatibility(sufficient) statistics" for model parameters.
                Default is \code{NULL} so the fitting will be based on
                the data. If the \code{CompStat} is provided then the inputs in
                \code{Data.Cur} and \code{Data.Hist} will be ignored.

     \code{n0} is the sample size of historical data.

     \code{mean0} is the sample mean of the historical data.

     \code{var0} is the sample variance of the historical data.

     \code{n1} is the sample size of current data.

     \code{mean1} is the sample mean of the current data.

     \code{var1} is the sample variance of the current data.
}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{(\mu, \sigma^2)} and \eqn{\delta}.
             The form of the prior for model parameter \eqn{(\mu, \sigma^2)} is \eqn{(1/\sigma^2)^a}.
             When \eqn{a = 1} it corresponds to the reference prior, and when \eqn{a = 1.5} it
             corresponds to the Jeffrey's prior.

     \code{joint.a} is the power \eqn{a} in formula \eqn{(1/\sigma^2)^a}, the prior for \eqn{(\mu, \sigma^2)} jointly.

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
\item{mu}{posterior of the model parameter \eqn{\mu}.}

\item{sigmasq}{posterior of the model parameter \eqn{\sigma^2}.}

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
NormalData0 <- rnorm(n = 100, mean= 20, sd = 1)

set.seed(12345)
NormalData1 <- rnorm(n = 50, mean= 30, sd = 1)

NormalNPP_MCMC(Data.Cur = NormalData1, Data.Hist = NormalData0,
               CompStat = list(n0 = 100, mean0 = 10, var0 = 1,
               n1 = 100, mean1 = 10, var1 = 1),
               prior = list(joint.a = 1.5, delta.alpha = 1, delta.beta = 1),
               MCMCmethod = 'RW', rw.logit.delta = 1,
               ind.delta.alpha= 1, ind.delta.beta= 1, nsample = 10000,
               control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
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
}

\keyword{normalized power prior}
