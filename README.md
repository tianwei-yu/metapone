The metapone package conducts pathway tests for untargeted metabolomics data. It has three main characteristics: (1) expanded database combining SMPDB and Mummichog databases, with manual cleaning to remove redundancies; (2) A new weighted testing scheme to address the issue of metabolite-feature matching uncertainties; (3) Can consider positive mode and negative mode data in a single analysis. 

It can be install by calling devtools::install_github("tianwei-yu/metapone"). To use the package, you need to have testing results on untargetted metabolomics data ready. The test result should contain at least three clumns - m/z, retention time, and feature p-value. An example input data can be seen here:

```{r}
> library(metapone)
> data(pos)
> head(pos)


       m.z retention.time    p.value statistic
1 85.04027       55.66454 0.22109229 -1.231240
2 85.07662       56.93586 0.52181695 -0.642790
3 85.57425      125.97117 0.13483680 -1.507372
4 86.06064      194.81306 0.26069118  1.131101
5 86.08001       54.74512 0.17206535  1.375352
6 86.09704      177.73650 0.07541608  1.796427
```

The package can use data from a single ion mode. But it can also handle the situation where data are collected on the same samples using both modes. If both positive mode and negative mode data are present, each is input into the algorithm as a separate matrix.

```{r example input second matrix}
> data(neg)
> head(neg)
               mz       chr      pval t-statistic
result.1 85.00448 268.83027 0.2423777   1.1690645
result.2 87.00881  48.84882 0.2222984   1.2204394
result.3 87.92531 161.99560 0.1341622   1.4978887
result.4 88.00399 129.88520 0.2941855  -1.0489839
result.5 88.01216  35.81698 0.8510984  -0.1877171
result.6 88.98808 127.47973 0.1748255  -1.3568608
```

The test is based on HMDB identification. The common adduct ions are pre-processed and stored in:

```{r example load database}
> data(hmdbCompMZ)
> head(hmdbCompMZ)
      HMDB_ID      m.z ion.type
1 HMDB0059597 1.343218     M+3H
2 HMDB0001362 1.679159     M+3H
3 HMDB0037238 2.341477     M+3H
4 HMDB0005949 3.345944     M+3H
5 HMDB0002387 4.011337     M+3H
6 HMDB0002386 4.677044     M+3H
```
Pathway information is built-in. It combines SMPDB and Mummichog databases. The two databases can be found here:

https://smpdb.ca/

https://shuzhao-li.github.io/mummichog.org/software.html

Here is the combined database after manual pruning of some highly-overlaping pathways:

```{r example load pathway}
> data(pa)
> head(pa)
      database         pathway.name     HMDB.ID KEGG.ID  category
1 Metapone 191 Pterine Biosynthesis HMDB0006822  C05922 Metabolic
2 Metapone 191 Pterine Biosynthesis HMDB0002111  C00001 Metabolic
3 Metapone 191 Pterine Biosynthesis HMDB0006821  C05923 Metabolic
4 Metapone 191 Pterine Biosynthesis HMDB0000142  C00058 Metabolic
5 Metapone 191 Pterine Biosynthesis HMDB0015532         Metabolic
6 Metapone 191 Pterine Biosynthesis HMDB0001273  C00044 Metabolic
```

The user can specify which adduct ions are allowed by setting the allowed adducts. For example:

```{r example adduct ions}
pos.adductlist = c("M+H", "M+NH4", "M+Na", "M+ACN+H", "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H")
neg.adductlist = c("M-H", "M-2H", "M-2H+Na", "M-2H+K", "M-2H+NH4", "M-H2O-H", "M-H+Cl", "M+Cl", "M+2Cl")
```

It is common for a feature to be matched to multiple metabolites. Assume a feature is matched to n metabolites, metapone weighs the feature by (1/n)^p, where p is a power term to tune the penalty. n can also be capped at a certain level such that the penalty is limited. These are controlled by parameters:

Setting p: fractional.count.power = 0.5

Setting the cap of n: max.match.count = 10 (this is to cap the level of penalty on a multiple-matched feature.

Other parameters include p.threshold, which controls which metabolic feature is considered significant. The testing is done by permutation. Overall, the analysis is conducted this way:

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
