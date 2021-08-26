# metapone
R package of pathway testing for untargeted metabolomics data
The metapone package conducts pathway tests for untargeted metabolomics data. It has three main characteristics: (1) expanded database combining SMPDB and Mummichog databases, with manual cleaning to remove redundancies; (2) A new weighted testing scheme to address the issue of metabolite-feature matching uncertainties; (3) Can consider positive mode and negative mode data in a single analysis. 

```{r setup}
library(metapone)
```

The input should contain at least three clumns - m/z, retention time, and feature p-value:

```{r example input}
data(pos)
head(pos)
```

If both positive mode and negative mode data are present, each is input into the algorithm as a separate matrix

```{r example input second matrix}
data(neg)
head(neg)
```

The test is based on HMDB identification. The common adduct ions are pre-processed and stored in:

```{r example load database}
data(hmdbCompMZ)
head(hmdbCompMZ)
```
Pathway information is built-in:

```{r example load pathway}
data(pa)
head(pa)
```

The user can specify which adduct ions are allowed by setting the allowed adducts. For example:

```{r example adduct ions}
pos.adductlist = c("M+H", "M+NH4", "M+Na", "M+ACN+H", "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H")
neg.adductlist = c("M-H", "M-2H", "M-2H+Na", "M-2H+K", "M-2H+NH4", "M-H2O-H", "M-H+Cl", "M+Cl", "M+2Cl")
```

It is common for a feature to be matched to multiple metabolites. Assume a feature is matched to n metabolites, metapone weighs the feature by (1/n)^p, where p is a power term to tune the penalty. n can also be capped at a certain level such that the penalty is limited. These are controlled by parameters:

Setting p: fractional.count.power = 0.5
Setting the cap of n: max.match.count = 10

Other parameters include p.threshold, which controls which metabolic feature is considered significant, and num.nodes, which controls the number of CPU cores to use in the computation. The testing is done by permutation. Overall, the analysis is conducted this way:

```{r example analysis}
r<-metapone(pos, neg, pa, hmdbCompMZ=hmdbCompMZ, pos.adductlist=pos.adductlist, neg.adductlist=neg.adductlist, p.threshold=0.05,n.permu=100,fractional.count.power=0.5, max.match.count=10)
hist(ptable(r)[,1])
```

We can subset the pathways that are significant:

```{r example continued}
selection<-which(ptable(r)[,1]<0.025)
ptable(r)[selection,]
ftable(r)[which(ptable(r)[,1]<0.025 & ptable(r)[,2]>=2)]
```
