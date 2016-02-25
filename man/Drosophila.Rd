\name{Drosophila}
\alias{Drosophila}
\docType{data}
\title{
Drosophila datasets and QTL mapping study
}
\description{
Two datasets with genotypes and phenotypes for backcrossed Drosophilas.
}
\usage{data(Drosophila)}
\format{
Two datasets with genotypes and phenotypes for backcrossed Drosophilas.

The set of 41 markers describes genotypes while 5 variables describe phenotypes. See references for more details.

\describe{
   \item{\code{bm}}{A data.frame with 370 obs. and 46 variables, first 41 are genotypes of gene markers, last five describes genotypes}
   \item{\code{bs}}{A data.frame with 402 obs. and 46 variables, first 41 are genotypes of gene markers, last five describes genotypes}
   \item{\code{chr}}{Factor w/ 4 levels CBFbeta, FLT3, None, Other}
   \item{\code{pos}}{Markers position on chromosom in centimorgnas}
 }
}
\details{
The phonotype pc1 is nicely described by genotype in both backcrossed datasets.
}
\source{
Genetic Architecture of a Morphological Shape Difference Between Two Drosophila Species
Zhao-Bang Zenga, Jianjun Liu, Lynn F. Stamb, Chen-Hung Kao, John M. Mercer, Cathy C. Laurie
Genetics, Vol. 154, 299-310, January 2000
}
\examples{
data(Drosophila)
library(lattice)
# calculate log likelihoods
pval1 = numeric(41)
for (i in 1:41) {
  y = Drosophila$bm$pc1
  x = factor(Drosophila$bm[,i])
  pval1[i] = logLik(lm(y~x))
}
# loglikelihood plot
xyplot(pval1~pos|chr, data=Drosophila, type=c("p","l"), 
     pch=19, ylab="log likelihood")
}
\keyword{Drosophila}

