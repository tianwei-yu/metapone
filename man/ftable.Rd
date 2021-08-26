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
\author{
  Tianwei Yu <yutianwei@cuhk.edu.cn>
}
\seealso{
  ptable
}
