\name{eunomia}
\alias{eunomia}
\docType{data}
\title{European Evaluation of Coercion in Psychiatry and Harmonisation of Best Clinical Practise}
\description{
This dataset touch one particular aspect from EUNOMIA dataset. The original dataset is much more richer.
}
\usage{data(eunomia)}
\format{
data.frame with 2008 obs. and 15 variables
\describe{
   \item{\code{CENTRE13}}{Center in which the patient is hospitalized, factor with 13 levels}
   \item{\code{SUBJECT}}{Patients ID}
   \item{\code{GENDER}, \code{AGE}, \code{NUM.HOSP}}{Gender, age and number of hospitalizations of given patient}
   \item{\code{CAT.T1}, \code{CAT.T2}, \code{CAT.T3}}{Clients Scale for Assessment of Treatment, short assessment, which measures the impact of COPD on a patients life, measured in times: T1, T2 and T3}
   \item{\code{BPRS.T1}, \code{BPRS.T2}, \code{BPRS.T3}}{Average score for Brief Psychiatric Rating Scale, measured in times: T1, T2 and T3}
   \item{\code{MANSA.T1}, \code{MANSA.T2}, \code{MANSA.T3}}{Scale which measures Quality of Life (Manchester Short Assessment of Quality of Life), measured in times: T1, T2 and T3}
   \item{\code{ICD10}}{International Statistical Classification of Diseases and Related Health Problems 10th Revision (ICD-10)}
  }
}
\details{
Artificial dataset generated to be consistent with Eunomia study (European Evaluation of Coercion in Psychiatry and Harmonisation of Best Clinical Practise) 
}
\source{
Artificial dataset generated to be consistent with Joanna R. study.

Eunomia dataset, http://www.eunomia-study.net/
}
\examples{
data(eunomia)
library(lattice)
bwplot(CENTRE13~BPRS.T1, data=eunomia)
xyplot(BPRS.T1~MANSA.T1|CENTRE13, data=eunomia, type=c("p","g","smooth"))
}
\keyword{eunomia}


