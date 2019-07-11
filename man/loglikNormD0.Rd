\name{loglikNormD0}
\alias{loglikNormD0}

\title{
A Function to Calculate Log-likelihood of the Historical Data, Given Array-valued Parameters, for Normal Population
}

\description{
The function returns a matrix of class "npp", each element is a log-likelihood of the historical data.
It is an intermediate step to calculate the "normalizing constant" \eqn{C(\delta)} in the normalized
power prior, for the purpose of providing a flexible implementation.
Users can specify their own likelihood function of the same class following this structure.
}


\usage{
loglikNormD0(D0, thetalist, ntheta = 2)
}

\arguments{
\item{D0}{a vector of each observation in historical data.}

\item{thetalist}{a list of parameter values. The number of elements is equal to \code{ntheta}.
                 Each element is a matrix. The sample should come from the posterior of the
                 powered likelihood for historical data, with each column corresponds to a
                 distinct value of the power parameter \eqn{\delta}
                 (the corresponding power parameter increases from left to right).
                 The number of rows is the number of Monte Carlo samples for each \eqn{\delta} fixed.
                 The number of columns is the number of selected knots (number of distinct \eqn{\delta}).
                 }

\item{ntheta}{a positive integer indicating number of parameters to be estimated in the model. }

}


\value{
A numeric matrix of log-likelihood, for the historical data given the matrix(or array)-valued parameters.
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
         \code{\link{loglikBerD0}};
         \code{\link{logCknot}};
         \code{\link{logCdelta}}
}

\keyword{normalized power prior}
