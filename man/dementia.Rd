\name{dementia}
\alias{dementia}
\docType{data}
\title{A set of datasets relatead to dementia}
\description{
Dataset from the book ,,Modele liniowe i mieszane w R, wraz z przykladami w analizie danych''.
}
\usage{data(dementia)}
\format{
data.frame with 1000 obs. and 4 variables
\describe{
   \item{\code{demscore}}{score of dementia}
   \item{\code{age}}{age, a factor with two levels}
   \item{\code{sex}}{sex, a factor with two levels}
   \item{\code{study}}{a source of data, a factor with 10 levels}
 }
}
\details{
Dataset from the book ,,Modele liniowe i mieszane w R, wraz z przykladami w analizie danych''.

Used as an example of mixed modeling in meta analysis.
}
\examples{
\dontrun{
  require(lme4)
  modelFullI <- lmer(demscore~age*sex+(age*sex|study), data=dementia,
                   REML=FALSE)
  summary(modelFullI)
}
}
\keyword{dementia}

