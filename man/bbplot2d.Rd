\name{bbplot2d}
\alias{bbplot2d}

\title{
Plot of metapone result.
}
\description{
The function \code{bbplot2d()} select important pathways with their P-value less than a threshold and returns a 2-D bubble plot with -log10(Pvalue) and the number of significant metabolites as coordinate axes.
}
\usage{
bbplot2d(res, p_thres = 0.05, sig_metab_thres = 1)
}

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
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ.metapone, p.threshold=0.05,
   n.permu=100,fractional.count.power=0.5, max.match.count=10)
bbplot2d(ptable(r)) # p_thres = 0.05
}
