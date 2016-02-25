\name{Plot results from post hoc testing}
\alias{plotPairwiseTests}
\title{
A function for visual representation of pairwise testing (both for pairwise.t.test and pairwise.wilcox.test)
}
\description{
Plot sets of groups in which means of medians are not significantly different.

On the veritical axis the means are marked. Then in a greedy fashion means that are not significantly different are linked by a line.
}
\usage{
 plotPairwiseTests(p.vals, means, alpha=0.05, digits=3, mar=c(2,10,3,1), ...) 
 
}

\arguments{
  \item{p.vals}{A slot \code{$p.value} from the object returned by pairwise.*.test function}
  \item{means}{A vector of means or medians corresponding to p.vals object (the order of groups should be the same in both objects)}
  \item{alpha}{A threshold for p.value}
  \item{digits}{Number of significant digits to be ploted with means.}
  \item{mar}{Figure margins, left margin should be large enought to handle names of groups}
  \item{\dots}{These arguments are passed to the plot function. }
}
\author{
Przemyslaw Biecek
}
\examples{
data(iris)
tmp1 <- pairwise.wilcox.test(iris$Sepal.Width, iris$Species)
tmp2 <- tapply(iris$Sepal.Width, iris$Species, median, na.rm=TRUE)
plotPairwiseTests(tmp1$p.value, tmp2, alpha=0.001)

}
