\name{ModeDeltaMultinomialNPP}
\alias{ModeDeltaMultinomialNPP}
\title{Calculate Posterior Mode of the Power Parameter in Normalized Power Prior with Grid Search, Multinomial Population}


\description{
The function returns the posterior mode of the power parameter \eqn{\delta} in multinomial population.
It calculates the log of the posterior density (up to a normalizing constant), and conduct a grid search
to find the approximate mode.
}


\usage{
ModeDeltaMultinomialNPP(Data.Cur, Data.Hist, CompStat = list(n0 = NULL, n1 = NULL),
                        npoints = 1000, prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5),
                        delta.alpha = 1, delta.beta = 1))
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

\item{npoints}{is a non-negative integer scalar indicating number of points on a regular spaced grid between [0, 1],
               where we calculate the log of the posterior and search for the mode.}

\item{prior}{a list of the hyperparameters in the prior for both \eqn{p} and \eqn{\delta}.

     \code{theta.dir} is a vector of \eqn{K} elements of the hyperparameter \eqn{\alpha}
                            in the prior distribution \eqn{Dir(\alpha[1],\alpha[2],...,\alpha[K])} for \eqn{\theta}.

     \code{delta.alpha} a scalar, the hyperparameter \eqn{\alpha} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.

     \code{delta.beta} a scalar, the hyperparameter \eqn{\beta} in the prior distribution \eqn{Beta(\alpha, \beta)} for \eqn{\delta}.
}
}


\value{
A numeric value between 0 and 1.
}


\details{
See example.
}



\examples{
ModeDeltaMultinomialNPP(CompStat = list(n0 = c(25,25,25,25), n1 = c(25,25,25,25)),
                        prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5, 0.5),
                                     delta.alpha = 1, delta.beta = 1))

ModeDeltaMultinomialNPP(CompStat = list(n0 = c(22,25,28,25), n1 = c(25,22,25,28)),
                        prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5, 0.5),
                                     delta.alpha = 1, delta.beta = 1))

ModeDeltaMultinomialNPP(CompStat = list(n0 = c(15,25,30,30), n1 = c(25,25,25,25)),
                        prior = list(theta.dir.alpha = c(0.5, 0.5, 0.5, 0.5),
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
         \code{\link{ModeDeltaPoisNPP}}
}

\keyword{normalized power prior}
