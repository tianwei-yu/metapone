\name{neg}
\alias{neg}
\docType{data}
\title{
Example negative mode data from the Metabolome Atlas of the Aging Mouse Brain
}
\description{
The data is generated from the hypocampus data of the Metabolome Atlas of the Aging Mouse Brain (ST001888) dataset. The p-values and test statistics were obtained by contrasting mouse hypocampus metabolome between prime-age mice and aging mice.}
\usage{data("neg")}
\format{
  A data frame with 6947 observations on the following 4 variables.
  \describe{
    \item{\code{m.z}}{a numeric vector. The mass-to-charge ratio of the features.}
    \item{\code{retention.time}}{a numeric vector. The retention time of the features.}
    \item{\code{p.value}}{a numeric vector. The p-values of the features.}
    \item{\code{statistic}}{a numeric vector. The test statistics of the features. }
	}
}
\references{
https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=FactorsData&StudyID=ST001888&StudyType=MS&ResultType=1
}
\examples{
data(neg)
}
\keyword{datasets}
