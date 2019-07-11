\name{VaccineData}
\alias{VaccineData}
\docType{data}
\title{
Dataset of a Vaccine Trial for RotaTeq and Multiple Historical Trials for Control Group
}
\description{
The study was designed to investigate the concomitant use of RotaTeq(Test Vaccine)
and some routine pediatric vaccines between 2001-2005.
The dataset includes four historical control trials.
The purpose of the study is to borrow the historical controls for the non-inferiority trial.
The interest is in the response rate to the routine vaccines.
}
\usage{data("VaccineData")}
\format{
  A data frame with 6 observations on the following 7 variables.
  \describe{
    \item{\code{Data.Time}}{indicator of historical data (coded as 0) or current data (coded as 1).}
    \item{\code{StudyID}}{character to distinguish different studies.}
    \item{\code{Group}}{indicator of control group (coded as 0) or treatment group (coded as 1). }
    \item{\code{Start.Year}}{start year of the trial }
    \item{\code{End.Year}}{end year of the trial}
    \item{\code{N}}{total number of patients enrolled and dosed in the group}
    \item{\code{y}}{total number of patients respond to the vaccine}
  }
}



\references{
Liu, G.F. (2018).
A Dynamic Power Prior for Borrowing Historical
Data in Noninferiority Trials with Binary Endpoint.
\emph{Pharmaceutical Statistics  17:61-73.}
}
\examples{
data(VaccineData)
}
\keyword{Bernoulli}
