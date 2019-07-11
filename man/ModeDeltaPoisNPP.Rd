\name{ModeDeltaPoisNPP}
\alias{ModeDeltaPoisNPP}
\title{Calculate Posterior Mode of the Power Parameter in Normalized Power Prior with Grid Search, Poisson Population}


\description{
The function returns the posterior mode of the power parameter \eqn{\delta} in multinomial population.
It calculates the log of the posterior density (up to a normalizing constant), and conduct a grid search
to find the approximate mode.
}


\usage{
ModeDeltaPoisNPP(Data.Cur, Data.Hist,
                 CompStat = list(n0 = NULL, mean0 = NULL, n1 = NULL, mean1 = NULL),
                 npoints = 1000, prior = list(lambda.shape = 1/2,
                 lambda.scale = 100, delta.alpha = 1, delta.beta = 1))
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

\item{npoints}{is a non-negative integer scalar indicating number of points on a regular spaced grid between [0, 1],
               where we calculate the log of the posterior and search for the mode.}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{\lambda} and \eqn{\delta}.
             A Gamma distribution is used as the prior of \eqn{\lambda}, and a Beta distribution is used as the prior of \eqn{\delta}.

     \code{lambda.shape} is the shape (hyper)parameter in the prior distribution \eqn{Gamma(shape, scale)} for \eqn{\lambda}.

     \code{lambda.scale} is the scale (hyper)parameter in the prior distribution \eqn{Gamma(shape, scale)} for \eqn{\lambda}.

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
ModeDeltaPoisNPP(CompStat = list(n0 = 50, mean0 = 10, n1 = 50, mean1 = 10), npoints = 1000,
                 prior = list(lambda.shape = 1/2, lambda.scale = 100,
                              delta.alpha = 1, delta.beta = 1))

ModeDeltaPoisNPP(CompStat = list(n0 = 50, mean0 = 10, n1 = 50, mean1 = 9.5), npoints = 1000,
                 prior = list(lambda.shape = 1/2, lambda.scale = 100,
                              delta.alpha = 1, delta.beta = 1))

ModeDeltaPoisNPP(CompStat = list(n0 = 50, mean0 = 10, n1 = 50, mean1 = 9), npoints = 1000,
                 prior = list(lambda.shape = 1/2, lambda.scale = 100,
                              delta.alpha = 1, delta.beta = 1))
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
         \code{\link{ModeDeltaBerNPP}};
         \code{\link{ModeDeltaNormalNPP}};
         \code{\link{ModeDeltaMultinomialNPP}}
}

\keyword{normalized power prior}
