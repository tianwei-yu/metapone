\name{metapone-package}
\alias{metapone-package}
\docType{package}
\title{
\packageTitle{metapone}
}
\description{
\packageDescription{metapone}
The package contains two types of statistical testing that considers matching uncertainty - (1) a permutation test that is based on the hypergeometric test and (2) a GSEA type test with weighted features/metabolites.

%The package conducts a permutation test that is based on the hypergeometric test for metabolic pathway analysis. 
}
\details{
%The package conducts a weighted hypergeometric test using permutations on metabolomics data. The weights are assigned based on how many metabolites each data feature can match to. 
The package conducts (1) a weighted hypergeometric test using permutations on metabolomics data. The weights are assigned based on how many metabolites each data feature can match to, (2) a GSEA type test based on an estimation of importance of metabolites/features. The importance is evluated by the size of matching for each metabolite/feature and the p-value of features.

The user can tune a parameter to change the penalty for multiple-matched features and choose the type of pathway testing. 
}
\author{
Tianwei Yu (\email{yutianwei@cuhk.edu.cn})
}
\keyword{ package }
