\name{SPDData}
\alias{SPDData}
\docType{data}
\title{
Dataset for Diagnostic Test (PartoSure Test, Medical Device) Evaluation for Spontaneous Preterm Delivery
}
\description{
The diagnostic test was developed to aid in rapidly assess the risk of
spontaneouspreterm delivery within 7 days from the time of diagnosis in pre-pregnant
women with signs and symptoms.
The same diagnostic test was used for two populations in US and EU respectively.
The number of counts in the four cells (True positive, false positive, false negative, true negative)
was recorded.
}
\usage{data("SPDData")}
\format{
  A data frame with 2 observations on the following 5 variables.
  \describe{
    \item{\code{Data.Region}}{region where the diagnostic test was conducted}
    \item{\code{TPDP}}{number of subjects with tested positive and the disease status positive (true positive)}
    \item{\code{TPDN}}{number of subjects with tested positive but the disease status negative (false positive)}
    \item{\code{TNDP}}{number of subjects with tested negative and the disease status positive (false negative)}
    \item{\code{TNDN}}{number of subjects with tested negative and the disease status negative (true negative)}
  }
}


\source{
https://www.accessdata.fda.gov/cdrh_docs/pdf16/P160052C.pdf
}

\examples{
data(SPDData)
}
\keyword{multinomial}
