\name{logCknot}
\alias{logCknot}

\title{
A Function to Calculate \eqn{logC(\delta)} on Selected Knots
}

\description{
The function returns a sequence of the values, \eqn{logC(\delta)} on selected knots, given
input vector of \eqn{\delta}.
}


\usage{
logCknot(deltaknot, llikf0)
}

\arguments{
\item{deltaknot}{a vector of the knots for \eqn{\delta}. It should be selected before conduct the sampling.}

\item{llikf0}{a matrix of the log-likelihoods of class "npp". }

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
         \code{\link{logCdelta}}
}

\keyword{normalized power prior}
