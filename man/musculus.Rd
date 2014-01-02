\name{musculus}
\alias{musculus}
\alias{muscle}
\docType{data}
\title{A dataset relatead to mice musculus growth which depends on diet and genetic structure}
\description{
Dataset from the book ,,Modele liniowe i mieszane w R, wraz z przykladami w analizie danych''.
}
\usage{data(musculus)}
\format{
data.frame with 30 obs. and 10 variables
\describe{
   \item{\code{id}}{an individual id}
   \item{\code{dadid}}{id of father, 0 for founders}
   \item{\code{momid}}{id of mother, 0 for founders}
   \item{\code{sex}}{sex}
   \item{\code{sigma}}{maximal stress}
   \item{\code{diet}}{diet, D1 or D2}
   \item{\code{k1}}{resilience coefficient in point 1}
   \item{\code{k2}}{resilience coefficient in point 2}
   \item{\code{E1}}{Younga module in point 1}
   \item{\code{E2}}{Younga module in point 2}
 }
}
\details{
Dataset from the book ,,Modele liniowe i mieszane w R, wraz z przykladami w analizie danych''.

Used as an example of model with mixed effects where random effects have know dependency structure, here related to the kinship coefficient.
}
\examples{
\dontrun{
   require(kinship2)
   pedmus <- pedigree(musculus$id, musculus$dadid, musculus$momid, musculus$sex)
   plot(pedmus, affected=musculus$diet)

   fam  <- makefamid(musculus$id, musculus$dadid, musculus$momid)
   kmatrix <- makekinship(fam, musculus$id, musculus$dadid, musculus$momid)
   kmatrix[1:5,1:15]
}
}
\keyword{musculus}

