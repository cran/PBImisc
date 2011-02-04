\name{kidney}
\alias{kidney}
\docType{data}
\title{Graft function after kidney transplantation}
\description{
Artificial dataset (subset from real dataset with some random modifications)
}
\usage{data(kidney)}
\format{
data.frame with 334 obs. and 16 variables
\describe{
   \item{\code{receiver.age}, \code{donor.age}}{Age od donor and receiver}
   \item{\code{CIT}}{Cold ischemia time}
   \item{\code{discrepancy.AB}, \code{discrepancy.DR}}{discrepancies in AB and DR antibodies}
   \item{\code{therapy}}{scheme of immunosuppression}
   \item{\code{diabetes}}{diabetes}
   \item{\code{bpl.drugs}}{number of drugs for blood pressure lowering}
   \item{\code{MDRD7}, \code{MDRD30}, \code{MDRD3}, \code{MDRD6}, \code{MDRD12}, \code{MDRD24}, \code{MDRD36}, \code{MDRD60}}{ MDRD (Modification of Diet in Renal Disease) as a estiamtor of glomerular filtration rate (GFR) from serum creatinine, measured 7, 30 days and 3, 6, 12, 24, 36 and 60 months after kidney transplantation}
 }
}
\details{
Example of longitudinal study, note that graft for all patients survives 5 years after kidney transplantation.
}
\source{
Artificial dataset generated to be consistent with Maria M. study
}
\examples{
data(kidney)
PBImisc:::boxplot.in.time(kidney[,9:16], colnames(kidney[,9:16]), additional=TRUE)
}
\keyword{kidney}

