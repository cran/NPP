\name{ModeDeltaLMNPP}
\alias{ModeDeltaLMNPP}
\title{Calculate Posterior Mode of the Power Parameter in Normalized Power Prior with Grid Search,
Normal Linear Model}


\description{
The function returns the posterior mode of the power parameter \eqn{\delta} in normal linear model.
It calculates the log of the posterior density (up to a normalizing constant), and conduct a grid search
to find the approximate mode.
}


\usage{
ModeDeltaLMNPP(y.Cur, y.Hist, x.Cur = NULL, x.Hist = NULL, npoints = 1000,
               prior = list(a = 1.5, b = 0, mu0 = 0, Rinv = matrix(1, nrow = 1),
                            delta.alpha = 1, delta.beta = 1))
}

\arguments{
\item{y.Cur}{a vector of individual level of the response y in current data.}

\item{y.Hist}{a vector of individual level of the response y in historical data.}


\item{x.Cur}{a vector or matrix or data frame of covariate observed in the current data.
             If more than 1 covariate available, the number of rows is equal to the number of
             observations.}

\item{x.Hist}{a vector or matrix or data frame of covariate observed in the historical data.
             If more than 1 covariate available, the number of rows is equal to the number of
             observations.}

\item{npoints}{is a non-negative integer scalar indicating number of points on a regular spaced grid between [0, 1],
               where we calculate the log of the posterior and search for the mode.}

\item{prior}{a list of the hyperparameters in the prior for model parameters
             \eqn{(\beta, \sigma^2)} and \eqn{\delta}.
             The form of the prior for model parameter \eqn{(\beta, \sigma^2)} is in the section "Details".

     \code{a} a positive hyperparameter for prior on model parameters. It is the power \eqn{a} in formula \eqn{(1/\sigma^2)^a};
     See details.

     \code{b} equals 0 if a flat prior is used for \eqn{\beta}. Equals 1 if a normal prior is used for \eqn{\beta}; See details.

     \code{mu0} a vector of the mean for prior \eqn{\beta|\sigma^2}. Only applicable if \code{b = 1}.

     \code{Rinv} inverse of the matrix \eqn{R}. The covariance matrix of the prior for \eqn{\beta|\sigma^2} is \eqn{\sigma^2 R^{-1}}.

     \code{delta.alpha} is the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

     \code{delta.beta} is the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.
}
}




\details{
If \eqn{b = 1}, prior for \eqn{(\beta, \sigma)} is \eqn{(1/\sigma^2)^a * N(mu0, \sigma^2 R^{-1})}, which includes the g-prior.
If \eqn{b = 0}, prior for \eqn{(\beta, \sigma)} is \eqn{(1/\sigma^2)^a}.
The outputs include posteriors of the model parameter(s) and power parameter, acceptance rate when sampling \eqn{\delta}, and
the deviance information criteria.
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
         \code{\link{ModeDeltaBerNPP}};
         \code{\link{ModeDeltaNormalNPP}};
         \code{\link{ModeDeltaMultinomialNPP}};
         \code{\link{ModeDeltaNormalNPP}}
}

\keyword{normalized power prior}
