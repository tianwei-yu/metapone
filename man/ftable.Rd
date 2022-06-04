\name{ftable}
\alias{ftable,metaponeResult-method}
\alias{ftable}
\title{
  Acessor functions for the feature mapping table in a
  metaponeResult object.
}
\description{
Returns a \code{list} containing the mapped features in each pathway.
}
\usage{
\S4method{ftable}{metaponeResult}(object)
}
\arguments{
  \item{object}{A metaponeResult object.}
}
\details{
Each pathway is represented by a \code{data.frame} as an item in the
		list object. The dataframe include information of m.z, 
		retention.time, p.value, statistic, HMDB_ID, theoretical m.z, 
		ion.type, fractional counts.
}
\value{
The method returns a list. Each item is for a pathway. Matched significant metabolites are included.
}
\author{
  Tianwei Yu <yutianwei@cuhk.edu.cn>
}
\seealso{
  ptable
}
\examples{
data(hmdbCompMZ.metapone)
data(pa)
data(pos)
data(neg)
dat <- list(pos, neg)
type <- list("pos", "neg")
r<-metapone(dat, type, pa, hmdbCompMZ=hmdbCompMZ.metapone, p.threshold=0.05,
   n.permu=100,fractional.count.power=0.5, max.match.count=10)
ftable(r)[1:6]
}

