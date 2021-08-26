\name{metaponeResult-class}
\docType{class}
\alias{metaponeResult}
\alias{metaponeResult-class}
\alias{metaponeResult-method}
\alias{ptable, metaponeResult-method}
\alias{ftable, metaponeResult-method}


\title{Class "metaponeResult"}
\description{	

  This class represents the results of pathway testing. The testing result 
  contain two major components: the significant level of each pathway, and 
  the features matched to each pathway.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("metaponeResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{test.result}:}{a dataframe containing p_value, n_significant 
	metabolites, n_mapped_metabolites, n_metabolites, significant metabolites, 
	mapped_metabolite IDs and pathway name.}
    \item{\code{mapped.features}:}{A list containing n entries, where n is the 
	number of pathways. Each entry is a data frame, containing the features 
	mapped to this pathway. The information include m.z, retention.time, 
	p.value, statistic, HMDB_ID, theoretical m.z, ion.type, fractional counts.}
  }
}

\section{Methods}{
  \describe{
    \item{ptable}{\code{signature(object = "metaponeResult")}:
        return the \code{data.frame} of test statistics for each pathway, 
		including p_value, n_significant metabolites, n_mapped_metabolites, 
		n_metabolites, significant metabolites, mapped_metabolite IDs 
		and pathway name.}

  \item{ftable}{\code{signature(object = "metaponeResult")}: Returns a
      \code{list} containing the mapped features in each pathway. Each
		pathway is represented by a \code{data.frame} as an item in the
		list object. The dataframe include information of m.z, 
		retention.time, p.value, statistic, HMDB_ID, theoretical m.z, 
		ion.type, fractional counts.}
  }
}
\author{Tianwei Yu}
\keyword{classes}
