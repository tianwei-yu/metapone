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
\author{
  Tianwei Yu <yutianwei@cuhk.edu.cn>
}
\seealso{
  ftable
}
