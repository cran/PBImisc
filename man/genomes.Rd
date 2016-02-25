\name{genomes}
\alias{genomes}
\docType{data}
\title{724 bacterial genomes data}
\description{
Few parameters gathered for 724 bacterial species.
}
\usage{data(genomes)}
\format{
data.frame with 724 obs. and 7 variables
\describe{
   \item{\code{organism}}{Organism name, unique value for every row}
   \item{\code{group}}{Group, a factor with 22 levels}
   \item{\code{size}}{Genome size in Mbp}
   \item{\code{CG}}{GC content for genome sequence}
   \item{\code{habitat}, \code{temp.group}, \code{temperature}}{Where does this bacteria live?}
  }
}
\details{
 This dataset is prepared by Pawel M., data are taken from NCBI repository.
 
 See http://www.ncbi.nlm.nih.gov/ for more details
}
\source{
Pawel M. study
}
\examples{
data(genomes)
library(ggplot2)
# is this relation linear ?
qplot(size,GC, data=genomes) + theme_bw()
# or linear in log scales?
qplot(size,GC, data=genomes, log="xy") + theme_bw()
}
\keyword{genomes}
