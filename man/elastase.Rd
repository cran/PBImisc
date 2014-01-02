\name{elastase}
\alias{elastase}
\docType{data}
\title{Relation between graft function and elastase}
\description{
Relation between graft function and elastase from nephrology study.
}
\usage{data(elastase)}
\format{
data.frame with 54 obs. and 5 variables
\describe{
   \item{\code{sex}, \code{age}, \code{weight}}{Patient's sex, age and weight}
   \item{\code{elastase}}{Elastase concentration}
   \item{\code{GFR}}{Patient's GFR (glomerular filtration rate)}
 }
}
\details{
Artificial dataset (real one with some random modifications). Do not use it for medical reasoning.
}
\source{
Artificial dataset generated to be consistent with Malgorzata L. study
}
\examples{
data(elastase)
library(lattice)
xyplot(GFR~elastase, data=elastase, type=c("p","r","g"))
}
\keyword{elastase}

