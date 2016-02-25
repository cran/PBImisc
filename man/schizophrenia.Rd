\name{schizophrenia}
\alias{schizophrenia}
\docType{data}
\title{Genetic backgroud of schizophrenia}
\description{
Dataset with genotypes and phenotypes for 98 patients with schizophrenia disorder. 
}
\usage{data(schizophrenia)}
\format{
data.frame with 98 obs. and 9 variables
\describe{
   \item{\code{NfkB}, \code{CD28}, \code{IFN}}{Genotypes for SNP mutations in selected three genes}
   \item{\code{Dikeos.manic}, \code{Dikeos.reality.distortion}, \code{Dikeos.depression}, \code{Dikeos.disorganization}, \code{Dikeos.negative}}{Dikeos scores for schizophrenia measured in five domains }
   \item{\code{Dikeos.sum}}{Sum of Dikeos scores}
}
}   
\details{
Alleles for two SNPs in genes: Nuclear Factor-Kappa Beta (NfkB) and Cluster of Differentiation 28 (CD28) were examined as well as mental health described by five scales (see Dikeos 2008 for more details).
}
\source{
Artificial dataset generated to be consistent with Dorota F. study
}
\examples{
data(schizophrenia)
attach(schizophrenia)
interaction.plot(CD28, NfkB, Dikeos.sum)
interaction.plot(NfkB, CD28, Dikeos.sum)
model.tables(aov(Dikeos.sum~NfkB*CD28))
}
\keyword{schizophrenia}

