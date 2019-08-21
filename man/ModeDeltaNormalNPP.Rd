\name{ModeDeltaNormalNPP}
\alias{ModeDeltaNormalNPP}
\title{Calculate Posterior Mode of the Power Parameter in Normalized Power Prior with Grid Search, Normal Population}


\description{
The function returns the posterior mode of the power parameter \eqn{\delta} in multinomial population.
It calculates the log of the posterior density (up to a normalizing constant), and conduct a grid search
to find the approximate mode.
}


\usage{
ModeDeltaNormalNPP(Data.Cur, Data.Hist,
                   CompStat = list(n0 = NULL, mean0 = NULL, var0 = NULL,
                                   n1 = NULL, mean1 = NULL, var1 = NULL),
                   npoints = 1000,
                   prior = list(a = 1.5, delta.alpha = 1, delta.beta = 1))
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

\item{npoints}{is a non-negative integer scalar indicating number of points on a regular spaced grid between [0, 1],
               where we calculate the log of the posterior and search for the mode.}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{(\mu, \sigma^2)} and \eqn{\delta}.
             The form of the prior for model parameter \eqn{(\mu, \sigma^2)} is \eqn{(1/\sigma^2)^a}.
             When \eqn{a = 1} it corresponds to the reference prior, and when \eqn{a = 1.5} it
             corresponds to the Jeffrey's prior.

     \code{a} is the power \eqn{a} in formula \eqn{(1/\sigma^2)^a}, the prior for \eqn{(\mu, \sigma^2)} jointly.

     \code{delta.alpha} is the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

     \code{delta.beta} is the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.
}
}


\value{
A numeric value between 0 and 1.
}




\details{
See example.
}



\examples{

ModeDeltaNormalNPP(CompStat = list(n0 = 50, mean0 = 0, var0 = 1,
                                   n1 = 50, mean1 = 0, var1 = 1), npoints = 1000,
                   prior = list(a = 1.5, delta.alpha = 1, delta.beta = 1))

ModeDeltaNormalNPP(CompStat = list(n0 = 50, mean0 = 0, var0 = 1,
                                   n1 = 40, mean1 = 0.2, var1 = 1), npoints = 1000,
                   prior = list(a = 1.5, delta.alpha = 1, delta.beta = 1))

ModeDeltaNormalNPP(CompStat = list(n0 = 50, mean0 = 0, var0 = 1,
                                   n1 = 40, mean1 = 0.6, var1 = 1), npoints = 1000,
                   prior = list(a = 1.5, delta.alpha = 1, delta.beta = 1))

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
         \code{\link{ModeDeltaMultinomialNPP}};
         \code{\link{ModeDeltaPoisNPP}}
}

\keyword{normalized power prior}
