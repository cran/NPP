\name{LaplacelogC}
\alias{LaplacelogC}

\title{
A Function to Calculate \eqn{logC(\delta)} Based on Laplace Approximation
}

\description{
The function assumes that the prior of the model parameters is very flat that had very minor impact on the
shape of the power prior (posterior based on the D0).
}

\usage{
LaplacelogC(delta, loglikmle, detHessian, ntheta)
}

\arguments{

\item{delta}{the power parameter between 0 and 1. The function returns \eqn{logC(\delta)} }

\item{loglikmle}{a scalar; the loglikelihood of the historical data evaluated at the maximum likelihood
estimates based on the historical data}

\item{detHessian}{determinant of the Hessian matrix evaluated at the loglikelihood function with respect to the maximum likelihood
estimates based on the historical data}

\item{ntheta}{an positive integer indicating number of parameters in the model }

}


\value{
\eqn{logC(\delta)} based on the Laplace approximation. Can be used for the posterior sampling in the normalized power prior.
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
         \code{\link{logCknot}}
}

\keyword{normalized power prior}
