\name{milk}
\alias{milk}
\docType{data}
\title{Milk yield data}
\description{
Milk yield data for 10 unrelated cows
}
\usage{data(milk)}
\format{
data.frame with 40 obs. and 2 variables
\describe{
   \item{\code{cow}}{cow id, a factor with 10 levels}
   \item{\code{milk.amount}}{milk amount in kgs per week}
 }
}
\details{
Weekly milk yield amount for 10 cows. For every cow 5 measurements are taken.
}
\examples{
data(milk)
library(lattice)
# change the order of levels
milk$cow = reorder(milk$cow, milk$milk.amount, mean)
#plot it
dotplot(cow~milk.amount, data=milk)
}
\keyword{milk}

