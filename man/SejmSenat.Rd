\name{SejmSenat}
\alias{SejmSenat}
\docType{data}
\title{
SejmSenat
}
\description{
Changes in word usage in consecutive Sejm and Senate cadencies
}
\usage{data(SejmSenat)}
\format{
contingency matrix with 973 27 rows and 8 columns
\describe{
   \item{\code{Sejm.I}, \code{Sejm.II}, \code{Sejm.III}, \code{Sejm.IV}, }{summary of records from four Sejm cadencies}
   \item{\code{Senat.II}, \code{Senat.III}, \code{Senat.IV}, \code{Senat.V}, }{summary of records from four Senate cadencies}
   \item{\code{adj},   \code{adja},  \code{adjp},  \code{adv},   \code{aglt},  \code{bedzie},\code{conj},  \code{depr},  \code{fin},   \code{ger},   \code{ign},   \code{imps},  \code{impt},  \code{inf},   \code{interp},\code{num},   \code{pact},  \code{pant},  \code{pcon},  \code{ppas},  \code{praet}, \code{pred},  \code{prep},  \code{qub},   \code{siebie},\code{subst}, \code{winien}}{word modes}
}
}
\details{
Word usage statistics generated from Sejm and Senat records
}
\source{
The IPI PAN Corpus
webpage http://korpus.pl/ 
}
\examples{
data(SejmSenat)
library(ca)
# can you see some patterns?
plot(ca(SejmSenat[-15,]), mass =c(TRUE,TRUE), arrows =c(FALSE,TRUE))
}
\keyword{SejmSenat}

