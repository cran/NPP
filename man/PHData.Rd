\name{PHData}
\alias{PHData}
\docType{data}
\title{
PH Data on four sites in Virginia
}
\description{
The dataset is used to assess if there is site impairment.
The site impairment is defined as whether the pH values at a site indicate that the site violates a
(lower) standard of 6.0 more than 10\% of the time.
}
\usage{data("PHData")}
\format{
  A data frame with 325 observations on the following 3 variables.
  \describe{
    \item{\code{Station}}{the site number, labeled as 1 to 4}
    \item{\code{Data.Time}}{indicator of historical data (coded as 0) or current data (coded as 1)}
    \item{\code{PH}}{value of PH on the site }
  }
}


\examples{
data(PHData)
}
\keyword{normal}
