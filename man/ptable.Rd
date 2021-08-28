\name{ptable}
\alias{ptable,metaponeResult-method}
\alias{ptable}
\title{
  Acessor functions for the test result table in a
  metaponeResult object.
}
\description{
return the \code{data.frame} of test statistics for each pathway.
}
\usage{
\S4method{ptable}{metaponeResult}(object)
}
\arguments{
  \item{object}{A metaponeResult object.}
}

\details{
		Includes p_value, n_significant metabolites, n_mapped_metabolites, 
		n_metabolites, significant metabolites, mapped_metabolite IDs 
		and pathway name. 
}
\value{
The method returns a data frame with 6 columns: "p_value", "n_significant metabolites", "n_mapped_metabolites", "n_metabolites", "significant metabolites", "mapped_metabolites".
 }
\author{
  Tianwei Yu <yutianwei@cuhk.edu.cn>
}
\seealso{
  ftable
}
\examples{
data(hmdbCompMZ.metapone)
data(pa)
data(pos)
data(neg)
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ.metapone, 
   p.threshold=0.05,n.permu=100,fractional.count.power=0.5, max.match.count=10)
head(ptable(r))
}