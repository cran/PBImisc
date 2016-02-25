\name{endometriosis}
\alias{endometriosis}
\docType{data}
\title{Endometriosis study}
\description{
How the endometriosis affects concetration of alpha and beta factors in the blood.
}
\usage{data(endometriosis)}
\format{
data.frame with 165 obs. and 4 variables
\describe{
   \item{\code{disease}}{disease, blood samples were taken from women with endometriosis of from healthy ones}
   \item{\code{phase}}{phase in the menstrual cycle as the examination day (proliferative or secretory)}
   \item{\code{alpha.factor}, \code{beta.factor}}{concentration of alpha and beta factors in blood}
  }
}
\details{
Dataset used as example of ANCOVA 
}
\source{
Artificial dataset generated to be consistent with Ula S. study
}
\examples{
data(endometriosis)
library(lattice)
xyplot(log(alpha.factor)~log(beta.factor)|disease*phase, 
            data=endometriosis, type=c("p", "r"))
summary(aov(alpha.factor~beta.factor*disease*phase, data=endometriosis))
}
\keyword{endometriosis}
