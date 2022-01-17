\name{pa}
\alias{pa}
\docType{data}
\title{
Pathway-metabolite match file.
}
\description{
mapps pathways with metabolites.
}
\usage{data("pa")}
\format{
  A data frame with 5395 observations on the following 5 variables.
  \describe{
    \item{\code{database}}{a character vector}
    \item{\code{pathway.name}}{a character vector}
    \item{\code{HMDB.ID}}{a character vector}
    \item{\code{KEGG.ID}}{a character vector}
    \item{\code{category}}{a character vector}
  }
}
\source{
\href{https://smpdb.ca/}{Small Molecule Pathway Database}

\href{https://shuzhao-li.github.io/mummichog.org/software.html}{Mummichog}
}
\examples{
data(pa)
}
\keyword{datasets}
