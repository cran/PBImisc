\name{boxplotpp}
\alias{boxplotpp}
\alias{boxplotInTime}
\title{boxplot plus plus }
\description{
boxplotpp
}
\usage{
boxplotpp(x, xname=seq(1:ncol(x)), utitle="", addLines=TRUE, 
           color = ifelse(addLines, "white","lightgrey"), ...) 

boxplotInTime(x, xname, additional=T, color = ifelse(additional, 
     "white","lightgrey"), main="", ylim=range(unlist(x),na.rm=T), ..., 
     points = dim(x)[2], at = 1:points)
}
\arguments{
  \item{x}{  TODO }
  \item{xname}{ TODO }
  \item{utitle}{ TODO }
  \item{addLines}{ TODO }
  \item{color}{ TODO }
  \item{additional}{ TODO }
  \item{main}{ TODO }
  \item{points}{ TODO }
  \item{at}{ TODO }
  \item{ylim}{ TODO }
  \item{...}{ TODO }
}
\details{
TODO
}
\value{
TODO
}
\author{ Przemyslaw Biecek}
\examples{
#TODO
}
\keyword{boxplotpp}
