\name{vaccination}
\alias{vaccination}
\docType{data}
\title{Effective dose study}
\description{
What is the minimal dose that is effective?
}
\usage{data(vaccination)}
\format{
data.frame with 100 obs. and 2 variables
\describe{
   \item{\code{response}}{a reaction effect}
   \item{\code{dose}}{a dose that was applied}
 }
}
\details{
Responses for different doses of treatment.
}
\source{
Artificial dataset generated to be consistent with Karolina P. study
}
\examples{
data(vaccination)
library(lattice)
bwplot(response~dose, data=vaccination)
}
\keyword{vaccination}

