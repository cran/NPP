\name{logCdelta}
\alias{logCdelta}

\title{
A Function to Interpolate \eqn{logC(\delta)} Based on Its Values on Selected Knots
}

\description{
The function returns the interpolated value (a scalar) of \eqn{logC(\delta)} based on its results on selected knots,
given input vector of \eqn{\delta}.
}


\usage{
logCdelta(delta, deltaknot, lCknot)
}

\arguments{

\item{delta}{a scalar of the input value of \eqn{\delta}.}

\item{deltaknot}{a vector of the knots for \eqn{\delta}. It should be selected before conduct the sampling.}

\item{lCknot}{a vector of the values \eqn{logC(\delta)} on selected knots, coming from the function \code{logCknot}.}

}


\value{
A sequence of the values, \eqn{logC(\delta)} on selected knots.
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
         \code{\link{loglikNormD0}};
         \code{\link{loglikBerD0}};
         \code{\link{logCknot}}
}

\keyword{normalized power prior}
