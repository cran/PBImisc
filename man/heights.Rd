\name{heights}
\alias{heights}
\docType{data}
\title{
Husband and Wife heights
}
\description{
A dataset from ,,A modern approach to regression with R''. Simon J. Sheather 2009 .
Paired heights for husbands and wifes.
}
\usage{data(heights)}
\format{
data.frame with 96 obs. and 2 variables
\describe{
   \item{\code{Husband}, \code{Wife}}{Height of husband and wife.}
}
}
\details{
The dataset from ,,A modern approach to regression with R''. Simon J. Sheather 2009 
}
\source{
A modern approach to regression with R. Simon J. Sheather 2009 
}
\examples{
data(heights)
plot(Husband~Wife, data=heights, pch=19)
abline(lm(Husband~Wife, data=heights), col="red")
abline(lm(Husband~Wife-1, data=heights), col="blue")
}
\keyword{heights}

