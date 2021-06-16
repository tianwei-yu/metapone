\name{hmdbCompMZ.metapone}
\alias{hmdbCompMZ.metapone}
\docType{data}
\title{
the m/z values of common adduct ions of metapone metaboites
}
\description{
Monoisotopic mass of common adduct ions, limited to those included in the pathways in metapone.
}
\usage{data("hmdbCompMZ.metapone")}
\format{
  A data frame with 79350 observations on the following 3 variables.
  \describe{
    \item{\code{HMDB_ID}}{HMDB ID.}
    \item{\code{ion.type}}{Adduct ion type.}
    \item{\code{m.z}}{the m/z of the adduct ion.}
  }
}
\details{
The main difference of using this dataset vs using hmdbCompMZ, is the metabolite universe in testing is limited to those metabolites matched to metapone pathways, not all HMDB metabolites.
}
\source{
https://hmdb.ca/
}
\references{
https://hmdb.ca/
}
\examples{
data(hmdbCompMZ)
}
\keyword{datasets}
