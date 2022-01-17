\name{neg}
\alias{neg}
\docType{data}
\title{
Example negative mode data from the Dorm Room Inhalation to Vehicle Emission (DRIVE) study
}
\description{
Data matrix from negative mode results. It is the test results from the data generated using C18 hydrophobic reversed-phase chromatography with negative ESI.
}
\usage{data("neg")}
\format{
  A data frame with 3337 observations on the following 4 variables.
  \describe{
    \item{\code{m.z}}{a numeric vector. The mass-to-charge ratio of the features.}
    \item{\code{retention.time}}{a numeric vector. The retention time of the features.}
    \item{\code{p.value}}{a numeric vector. The p-values of the features.}
    \item{\code{statistic}}{a numeric vector. The test statistics of the features. }
	}
}
\references{
Liang D, Moutinho JL, Golan R, Yu T, Ladva CN, Niedzwiecki M, Walker DI, Sarnat SE, Chang HH, Greenwald R, Jones DP, Russell AG, Sarnat JA (2018) Use of high-resolution metabolomics for the identification of metabolic T signals associated with traffic-related air pollution. Environment International. 120:145-154.
}
\examples{
data(neg)
}
\keyword{datasets}
