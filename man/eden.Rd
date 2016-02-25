\name{eden}
\alias{eden}
\docType{data}
\title{European day hospital evaluation}
\description{
This dataset bases on origical study of European day hospital evaluation

Artificial dataset (subset from real dataset with some random modifications). Do not use it for derivation of real conclusions.
}
\usage{data(eden)}
\format{
data.frame with 642 obs. and 12 variables
\describe{
   \item{\code{mdid}}{Medical doctor id, there are 24 different MDs which examine patients}
   \item{\code{center}}{City in which the examination takes place}
   \item{\code{BPRS.Maniac}, \code{BPRS.Negative}, \code{BPRS.Positive}, \code{BPRS.Depression}}{BPRS stands for Brief Psychiatric Rating Scale, scores are averaged in four subscales}
   \item{\code{BPRS.Average}}{Average from 24 questions}
   \item{\code{MANSA}}{Scale which measures Quality of Life (Manchester Short Assessment of Quality of Life)}
   \item{\code{sex}}{Sex}
   \item{\code{children}}{Number of childs}
   \item{\code{years.of.education}}{Number of years of education}
   \item{\code{day}}{Hospitalization mode, day or stationary}
 }
}
\details{
This dataset touch one particular aspect from EDEN dataset. The original dataset is much more richer.
}
\source{
Artificial dataset generated to be consistent with Joanna R. study.

Bases on European day hospital evaluation, http://www.edenstudy.com/
}
\examples{
data(eden)
library(lattice)
xyplot(BPRS.Average~MANSA|center, data=eden, type=c("p","g","smooth"))
}
\keyword{eden}
