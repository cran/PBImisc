\name{milkgene}
\alias{milkgene}
\docType{data}
\title{Mutation in BTN3A1 gene and milk yield}
\description{
It is known that BTN3A1 (Butyrophilin subfamily 3 member A1) has a crucial function in the secretion of lipids into milk. 
Doeas the SNP mutation in it change the average milk yield?

}
\usage{data(milkgene)}
\format{
data.frame with 1000 obs. and 5 variables
\describe{
   \item{\code{cow.id}}{cow id, there is 465 cows in this study}
   \item{\code{btn3a1}}{btn3a1 genotype, a factor with two levels}
   \item{\code{lactation}}{for some cows there are milk yileds for four lactations for other only for the first one}
   \item{\code{milk}, \code{fat}}{milk and fat amount in kgs per lactation}
 }
}
\details{
Milk and fat yields for 465 cows. For every cow also the genotype of btn3a1 is measured.
}
\source{
Artificial dataset generated to be consistent with Joanna Sz. study
}
\examples{
data(milkgene)
library(lattice)
xyplot(milk~fat, data=milkgene)
bwplot(milk~lactation, data=milkgene)
}
\keyword{milk}

