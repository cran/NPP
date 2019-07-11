\name{MultinomialNPP_MCMC}
\alias{MultinomialNPP_MCMC}
\title{MCMC Sampling for Multinomial Population using Normalized Power Prior}


\description{
Conduct posterior sampling for multinomial population with normalized power prior.
For the power parameter \eqn{\delta}, a Metropolis-Hastings algorithm with either
independence proposal, or a random walk proposal on its logit scale is used.
For the model parameter vector \eqn{\theta}, Gibbs sampling is used.
Assume the prior for model parameter \eqn{\theta} comes from a Dirichlet distribution.
}


\usage{
MultinomialNPP_MCMC(Data.Cur = c(10, 10, 10), Data.Hist = c(10, 10, 10),
                    CompStat = list(n0 = NULL, n1 = NULL),
                    prior = list(theta.dir = c(0.5, 0.5, 0.5),
                                 delta.alpha = 1, delta.beta = 1),
                    MCMCmethod = 'IND', rw.logit.delta = 0.1,
                    ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 5000,
                    control.mcmc = list(delta.ini = NULL, burnin = 0, thin = 1))
}

\arguments{
\item{Data.Cur}{a non-negative integer vector of \eqn{K} elements:
               c(number of success in group 1, number of success in group 2, ..., number of success in group \eqn{K})
               in the current data.}

\item{Data.Hist}{a non-negative integer vector of \eqn{K} elements:
               c(number of success in group 1, number of success in group 2, ..., number of success in group \eqn{K})
               in the historical data.}

\item{CompStat}{a list of two elements that represents the
                "compatibility(sufficient) statistics" for \eqn{\theta}.
                Default is \code{NULL} so the fitting will be based on
                the data. If the \code{CompStat} is provided then the inputs in
                \code{Data.Cur} and \code{Data.Hist} will be ignored.
                Note: in multinomial case providing \code{CompStat} is equivalent to provide the data summary as in
                \code{Data.Cur} and \code{Data.Cur}.

     \code{n0} is a non-negative integer vector of \eqn{K} elements for compatible statistics in historical data:
               c(number of success in group 1, number of success in group 2, ..., number of success in group \eqn{K}).

     \code{n1} is a non-negative integer vector of \eqn{K} elements for compatible statistics in current data:
               c(number of success in group 1, number of success in group 2, ..., number of success in group \eqn{K}).

}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{p} and \eqn{\delta}.

     \code{theta.dir} is a vector of \eqn{K} elements of the hyperparameter \eqn{\alpha}
                            in the prior distribution
                            \eqn{Dir(\alpha[1],\alpha[2],...,\alpha[K])} for \eqn{\theta}.

     \code{delta.alpha} a scalar, the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

     \code{delta.beta} a scalar, the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

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
\item{p}{posterior of the model parameter \eqn{\theta}.}

\item{delta}{posterior of the power parameter \eqn{\delta}.}

\item{acceptance}{the acceptance rate in MCMC sampling for \eqn{\delta} using Metropolis-Hastings algorithm. }

\item{DIC}{the deviance information criteria for model diagnostics.}

}





\details{
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate in sampling \eqn{\delta}, and
the deviance information criteria.
}



\examples{
MultinomialNPP_MCMC(Data.Cur = c(3,11,3,669), Data.Hist = c(9,20,9,473),
                    prior = list(theta.dir = c(1,1,1,1),
                    delta.alpha = 1, delta.beta = 1),
                    MCMCmethod = 'IND', rw.logit.delta = 0.1,
                    ind.delta.alpha = 1, ind.delta.beta = 1, nsample = 10000,
                    control.mcmc = list(delta.ini = NULL,
                    burnin = 2000, thin = 5))
}


\author{
Tianyu Bai \email{tianyu.bai24@gmail.com}
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
         \code{\link{BerNPP_MCMC}};
         \code{\link{NormalNPP_MCMC}};
         \code{\link{PoissonNPP_MCMC}}
}

\keyword{normalized power prior}
