\name{ecap}
\alias{ecap}
\docType{data}
\title{
Epidemiology of Allergic Disease in Poland 
}
\description{
This dataset touch one particular aspect from ECAP dataset. The original dataset is much more richer.
}
\usage{data(ecap)}
\format{
data.frame with 2102 obs. and 9 variables
\describe{
   \item{\code{city}, \code{district}}{City and district, city is a factor with nine levels, the district effect is nested in the city effect}
   \item{\code{sex}}{Sex}
   \item{\code{weight}, \code{height}}{Weight and height}
   \item{\code{house.surface}}{Surface of house in which the pearson live}
   \item{\code{PNIF}}{Peak Nasal Inspiratory Flow}
   \item{\code{age}}{Age of the pearson}
   \item{\code{allergenes}}{Number of allergens}
 }
}
\details{
PNIF stands for Peak Nasal Inspiratory Flow
}
\source{
Artificial dataset generated to be consistent with ECAP (Epidemiologia Chorob Alergicznych w Polsce) study
http://www.ecap.pl/
}
\examples{
data(ecap)
library(lattice)
xyplot(PNIF~age|city, data=ecap, type=c("p","g","smooth"))
}
\keyword{ecap}
