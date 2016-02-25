\name{apartments}
\alias{apartments}
\docType{data}
\title{Apartment prices in Warsaw in years 2007-2009}
\description{
Dataset downloaded from website http://www.oferty.net/. Dataset contains offer and transictional prices for apartments sold in in Warsaw in years 2007-2009.
}
\usage{data(apartments)}
\format{
data.frame with 973 obs. and 16 variables
\describe{
   \item{\code{year}}{data year of the transaction}
   \item{\code{month}}{data month of the transaction}
   \item{\code{surface}}{apartment area in m2}
   \item{\code{city}}{city (all transactions are from Warsaw)}
   \item{\code{district}}{district in which the apartment is located, factor with 28 levels}
   \item{\code{street}}{steet in which the apartment is located}
   \item{\code{n.rooms}}{number of rooms}
   \item{\code{floor}}{floor}
   \item{\code{construction.date}}{the construction year}
   \item{\code{type}}{ownership rights}
   \item{\code{offer.price}}{price in the offer}
   \item{\code{transaction.price}}{declared price in the transaction}
   \item{\code{m2.price}}{price per m2}
   \item{\code{condition}}{apartment condition, factor with 5 levels}
   \item{\code{lat}, \code{lon}}{latitude and longitude coordinates for district center}
 }
}
\details{
This and other related dataset you may find here http://www.oferty.net/. 
}
\source{
website http://www.oferty.net/
}
\examples{
data(apartments)
library(lattice)
xyplot(m2.price~construction.date|district, apartments, type=c("g","p"))

# 
# apartments2 = na.omit(apartments[,c(13,1,3,5,7,8,9,10,14,15,16)])
# wsp = (bincombinations(10)==1)[-1,]
# params = matrix(0, nrow(wsp), 3)
# for (i in 1:nrow(wsp)) {
# 	  model = lm(m2.price~., data=apartments2[,c(TRUE,wsp[i,])])
#   	params[i,1] = AIC(model, k=log(nrow(apartments2)))
#  	  params[i,2] = model$rank
#  	  params[i,3] = summary(model)$adj.r.squared
# }
# plot(params[,2], params[,3], xlab="no. of regressors", ylab="adj R^2")
# 
}
\keyword{apartments}

