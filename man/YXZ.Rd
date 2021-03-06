\name{YXZ}
\alias{YXZ}
\docType{data}
\title{
Artificial dataset which shows the differences between tests type I and III (sequential vs. marginal)
}
\description{
Artificial dataset, shows inconsistency for test type I and III
}
\usage{data(YXZ)}
\format{
data.frame with 100 obs. and 3 variables
\describe{
   \item{\code{X}, \code{Z}}{explanatory variables}
   \item{\code{Y}}{response variable}
}
}
\details{
See the example, results for staistical tests are inconsistet due to correlation between X and Z variables
}
\source{
Artificial dataset, generated by PBI
}
\examples{
attach(YXZ)
summary(lm(Y~X+Z))
anova(lm(Y~Z+X))
anova(lm(Y~X))
anova(lm(Y~Z))
}
\keyword{YXZ}

