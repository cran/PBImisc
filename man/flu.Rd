\name{flu}
\alias{flu}
\docType{data}
\title{
Numbers of flu occurences in the 10 years period in the Poland.
}
\description{
Data from National Institute of Hygiene reports.
Each row correspond to one record from NIH institute. 
}
\usage{data(flu)}
\format{
data.frame with 6384 obs. and 11 variables
\describe{
   \item{\code{region}}{Region for which given report was taken. A factor with 16 levels}
   \item{\code{inception.no}}{Number of flu occurences in given region for given report period (one or two weeks)}
   \item{\code{inception.no}}{Number of flu occurences in given region for given report period (one or two weeks)}
   \item{\code{inception.rate}}{Number of flu occurences normalized to 100k people}
   \item{\code{inception.no.0-14}, \code{inception.no.15+}, \code{inception.rate.0-14}, \code{inception.rate.15+}}{Absolute and normalized numbers of flu occurences calculated for age group 0-14 or 15+}
   \item{\code{date}}{Date of given report}
   \item{\code{date.id}}{Report id, there is 38 reports per year}
   \item{\code{latitude}, \code{longitude}}{Geographical coordinates for region}
 }
}
\details{
Dataset used during the third edition of WZUR conference, see http://www.biecek.pl/WZUR3/wzurDane.html for more information.
}
\source{
Reports from National Institute of Public Health - National Institute of Hygiene, see: http://www.pzh.gov.pl

More information: http://www.biecek.pl/WZUR3/wzurDane.html
}
\examples{
data(flu)
library(ggplot2)
subflu = flu[flu$region=="Mazowieckie", ]
# linear scale
qplot(date, inception.rate,data=subflu, geom="line")+scale_y_sqrt() +theme_bw()

# polar coordinates
qplot(1 + date.id*12/38, inception.rate,data=subflu, geom="path", xlab="month")+
        scale_y_sqrt()+geom_smooth(span=0.1,se=FALSE, size=2,col="red") + 
        coord_polar() +theme_bw()
}
\keyword{flu}

