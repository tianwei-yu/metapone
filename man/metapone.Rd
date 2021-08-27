\name{metapone}
\alias{metapone}
\title{
METAbolic pathway testing using both POsitive and NEgative mode data
}
\description{
Metapone conducts pathway tests for untargeted metabolomics data. It has three main characteristics: (1) expanded database combining SMPDB and Mummichog databases, with manual cleaning to remove redundancies; (2) A new weighted testing scheme to address the issue of metabolite-feature matching uncertainties; (3) Can consider positive mode and negative mode data in a single analysis.
}
\usage{
metapone(pos = NULL, neg = NULL, pa, hmdbCompMZ, pos.adductlist = 
c("M+H", "M+NH4", "M+Na", "M+ACN+H", "M+ACN+Na", "M+2ACN+H", "2M+H", 
"2M+Na", "2M+ACN+H"), neg.adductlist = c("M-H", "M-2H", "M-2H+Na", 
"M-2H+K", "M-2H+NH4", "M-H2O-H", "M-H+Cl", "M+Cl", "M+2Cl"), 
use.fractional.count = TRUE, match.tol.ppm = 5, p.threshold = 0.05, 
n.permu = 200, fractional.count.power = 0.5, max.match.count = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{
The positive ion mode test results. A matrix with four columns: m/z, retention time, p-value, test statistic. The package doesn't require both pos and neg to be present. One ion mode result is sufficient.
}
  \item{neg}{
The negative ion mode test results. A matrix with four columns: m/z, retention time, p-value, test statistic. The package doesn't require both pos and neg to be present. One ion mode result is sufficient.
}
  \item{pa}{
Pathway information. A data frame with five columns: database pathway ID, pathway name, HMDB ID, KEGG ID, category of pathway.
}
  \item{hmdbCompMZ}{
the m/z values of common adduct ions of HMDB metaboites. See the help file of hmdbCompMZ for details.
}
  \item{pos.adductlist}{
The vector of positive adduct ions to be considered.
}
  \item{neg.adductlist}{
The vector of negative adduct ions to be considered.
}
  \item{use.fractional.count}{
A lot of features match to multiple metabolites by m/z. Whether to discount such matches by using fractional counts.
}
  \item{match.tol.ppm}{
The ppm level when conducting m/z match.
}
  \item{p.threshold}{
The threshold of p-values of metabolic features to be considered significant.
}
  \item{n.permu}{
The number of permutations in permutation test.
}
  \item{fractional.count.power}{
The fractional counts are taken to this power to transform the weights.
}
  \item{max.match.count}{
When calculating fractional counts, some features might be matched to too many. In that case the number of matches is capped by the value of max.match.count.} 
}
\value{
The method returns a generic S4 object of class "metapone.result":
  \item{@test.results}{A matrix with 6 columns: "p_value", "n_significant metabolites", "n_mapped_metabolites", "n_metabolites", "significant metabolites", "mapped_metabolites". Each row is for a pathway.}
  \item{@mapped.features}{A list. Each item is for a pathway. The item lists matched significant metabolites.}
}
\references{
https://smpdb.ca/
https://shuzhao-li.github.io/mummichog.org/software.html
}
\author{
Tianwei Yu<yutianwei@cuhk.edu.cn>)
}
\seealso{
hmdbCompMZ, pa
}
\examples{
data(hmdbCompMZ.metapone)
data(pa)
data(pos)
data(neg)
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ.metapone, p.threshold=0.05,
   n.permu=100,fractional.count.power=0.5, max.match.count=10)
hist(ptable(r)[,1])
}