\name{bbplot1d}
\alias{bbplot1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of metapone result.
}
\description{
The function \code{bbplot1d()} select important pathways with their P-value less than a threshold and returns ranked bubble plot showing important pathways names and their corresponding -log10(Pvalue).
}
\usage{
bbplot1d(res, p_thres = 0.05, sig_metab_thres = 1, low.color = "MidnightBlue", high.color = "LightSkyBlue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{res}{
The result matrix obtained from metapone with columns: "p_value", "n_significant metabolites", "n_mapped_metabolites", "n_metabolites", "significant metabolites", "mapped_metabolites", "fdr".
}
\item{p_thres}{
The threshold of P-value for pathways to be shown in the bubble plot. The default threshold is 0.05.
}
\item{sig_metab_thres}{
The threshold of fractional matched significant metabolite count for pathways to be shown in the bubble plot. The default is 1. 
} 
\item{low.color}{
The GRB color of the lowest ldfr value to be shown in the bubble plot. 
}
\item{high.color}{
The GRB color of the highest ldfr value to be shown in the bubble plot. 
}
}
\author{
Leqi Tian (\email{leqitian@link.cuhk.edu.cn})
}
\seealso{
\code{\link{metapone}}
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
bbplot1d(ptable(r)) # p_thres = 0.05
}
