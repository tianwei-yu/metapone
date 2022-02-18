\name{metapone}
\alias{metapone}
\title{
METAbolic pathway testing using both POsitive and NEgative mode data
}
\description{
Metapone conducts pathway tests for untargeted metabolomics data. It has three main characteristics: (1) expanded database combining SMPDB and Mummichog databases, with manual cleaning to remove redundancies; (2) A new weighted testing scheme to address the issue of metabolite-feature matching uncertainties; (3) Can consider positive mode and negative mode data in a single analysis.
}
\usage{
metapone(pos=NULL, neg=NULL, pa, hmdbCompMZ, pos.adductlist = c("M+H", 
"M+NH4", "M+Na", "M+ACN+H","M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", 
"2M+ACN+H"), neg.adductlist = c("M-H","M-2H","M-2H+Na","M-2H+K", 
"M-2H+NH4","M-H2O-H","M-H+Cl", "M+Cl", "M+2Cl"),
use.fractional.count=TRUE, match.tol.ppm=5, p.threshold=0.05, 
n.permu=200, fractional.count.power=0.5, max.match.count=10, 
use.fgsea = FALSE, use.meta = FALSE)
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
  \item{use.fgsea}{
Whether to use a GSEA type test when performing pathway testing. When it is FALSE, a permutation-based weighted hypergeometric test is performed.
} 
  \item{use.meta}{
Whether to perform a GSEA type test with weighted metabolites. When it is FALSE, a GSEA type test is performed on weighted features.
}
}
\value{
The method returns a generic S4 object of class "metapone.result":
  \item{@test.results}{A matrix with 8 columns: "p_value", "n_significant metabolites", "n_mapped_metabolites", "n_metabolites", "significant metabolites", "mapped_metabolites", "lfdr", "adjust.p". Each row is for a pathway. When using GSEA test, "ES", "NES", "nMoreExtreme" are returned additionally.}
  \item{@mapped.features}{A list. Each item is for a pathway. The item lists matched significant metabolites.}
The columns in test.result are the following:
  \item{p_value}{
The p-value for each enrichment.
}
  \item{n_significant metabolites}{
The number of weighted significant metabolites associated with the enrichment.
}
  \item{n_mapped_metabolites}{
The number of weighted metabolites associated with the enrichment.
}
  \item{n_metabolites}{
The number of metabolites associated with the enrichment.  
}
  \item{significant metabolites}{
A string with the names of significant metabolites that drive the enrichment.
}
  \item{mapped_metabolites}{
A string with the names of metabolites that drive the enrichment.
}
  \item{lfdr}{
The local fdr value for each enrichment.
}
  \item{adjust.p}{
The enrichment BH-adjusted p-value for each enrichment.
}
  \item{ES}{
The enrichment score (Avaliable in GSEA test).
}
  \item{NES}{
The enrichment score normalized to mean enrichment of random samples of the same size (Avaliable in GSEA test).
}
  \item{nMoreExtreme}{
The number of times a random metabolite set had a more extreme enrichment score value (Avaliable in GSEA test).
}
}
\references{
\href{https://smpdb.ca/}{Small Molecule Pathway Database}

\href{https://shuzhao-li.github.io/mummichog.org/software.html}{Mummichog}
}
\author{
Tianwei Yu (\email{yutianwei@cuhk.edu.cn})
Leqi Tian (\email{leqitian@link.cuhk.edu.cn})
}
\seealso{
\code{\link{pa}}, \code{\link{hmdbCompMZ}}
}
\examples{
data(hmdbCompMZ.metapone)
data(pa)
data(pos)
data(neg)
# Permutation-based weighted hypergeometric test
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ.metapone, p.threshold=0.05,
   n.permu=100,fractional.count.power=0.5, max.match.count=10)
hist(ptable(r)[,1])

# Metabolites based GSEA test
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ.metapone, p.threshold=0.05,
   n.permu=100,fractional.count.power=0.5, max.match.count=10, use.fgsea = TRUE, use.meta = TRUE)
hist(ptable(r)[,1])

# Features based GSEA test
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ.metapone, p.threshold=0.05,
   n.permu=100,fractional.count.power=0.5, max.match.count=10, use.fgsea = FALSE, use.meta = FALSE)
hist(ptable(r)[,1])
}