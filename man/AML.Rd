\name{AML}
\alias{AML}
\docType{data}
\title{
Acute myeloid leukemia AML study
}
\description{
This dataset bases on blood samples for patients with Acute myeloid leukemia.
}
\usage{data(AML)}
\format{
data.frame with 66 obs. and 5 variables
\describe{
   \item{\code{Mutation}}{Factor w/ 4 levels CBFbeta, FLT3, None, Other}
   \item{\code{CD14.control}}{CD14 level in the control group}
   \item{\code{CD14.D3}}{CD14 level after D3 treatment}
   \item{\code{CD14.1906}}{CD14 level after D3 homolog 1906 treatment}
   \item{\code{CD14.2191}}{CD14 level after D3 homolog 2191 treatment}
 }
}
\details{
 Mutation - mutated gene that causes leucemia, one of following CBFbeta, FLT3, None, Other
 CD14.control, CD14.D3, CD14.1906, CD14.2191 - effects in vitamin D3 or its homologues
}
\source{
Artificial dataset generated to be consistent with Ewa M. study
}
\examples{
library(lattice)
data(AML)
AML2 = reshape(AML, direction="long", varying=colnames(AML)[2:5])
bwplot(CD14~time|Mutation, AML2)
interaction.plot(AML2$time,AML2$Mutation, AML2$CD14)
}
\keyword{AML}
