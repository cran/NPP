\name{ModeDeltaBerNPP}
\alias{ModeDeltaBerNPP}
\title{Calculate Posterior Mode of the Power Parameter in Normalized Power Prior with Grid Search, Bernoulli Population}


\description{
The function returns the posterior mode of the power parameter \eqn{\delta} in Bernoulli population.
It calculates the log of the posterior density (up to a normalizing constant), and conduct a grid search
to find the approximate mode.
}


\usage{
ModeDeltaBerNPP(Data.Cur, Data.Hist,
                CompStat = list(n0 = NULL, y0 = NULL, n1 = NULL, y1 = NULL),
                npoints = 1000,
                prior = list(p.alpha = 1, p.beta = 1,
                             delta.alpha = 1, delta.beta = 1))
}
\arguments{
\item{Data.Cur}{a non-negative integer vector of two elements: c(number of success, number of failure) in the current data.}

\item{Data.Hist}{a non-negative integer vector of two elements: c(number of success, number of failure) in the historical data.}

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

\item{npoints}{is a non-negative integer scalar indicating number of points on a regular spaced grid between [0, 1],
               where we calculate the log of the posterior and search for the mode.}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{p} and \eqn{\delta}.

     \code{p.alpha} is the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{p}.

     \code{p.beta} is the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{p}.

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
ModeDeltaBerNPP(Data.Cur = c(100, 40), Data.Hist = c(100, 40), npoints = 1000,
                prior = list(p.alpha = 1, p.beta = 1, delta.alpha = 1, delta.beta = 1))

ModeDeltaBerNPP(Data.Cur = c(100, 40), Data.Hist = c(100, 35), npoints = 1000,
                prior = list(p.alpha = 1, p.beta = 1, delta.alpha = 1, delta.beta = 1))

ModeDeltaBerNPP(Data.Cur = c(100, 40), Data.Hist = c(100, 50), npoints = 1000,
                prior = list(p.alpha = 1, p.beta = 1, delta.alpha = 1, delta.beta = 1))
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
         \code{\link{ModeDeltaNormalNPP}};
         \code{\link{ModeDeltaPoisNPP}};
         \code{\link{ModeDeltaMultinomialNPP}}
}

\keyword{normalized power prior}
