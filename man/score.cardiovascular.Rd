\name{score.cardiovascular}
\alias{score.cardiovascular}
\alias{calculateScoreEur}
\title{
SCORE for Cardiovascular Risk
}
\description{
Calculation of risk SCORE for use in the clinical management of cardiovascular risk in European.
}
\usage{
calculateScoreEur(age, cholesterol, SBP, currentSmoker, 
  gender = "Men", risk = "Low risk")

}
\arguments{
  \item{age}{  age in years }
  \item{cholesterol}{ in mmol/L }
  \item{SBP}{ Systolic blood pressure in mmHg }
  \item{currentSmoker}{ the current smoker status,  1 for current smokers, 0 for non smokers}
  \item{gender}{ "Men" or "Women" }
  \item{risk}{ is it "Low risk" or "High risk" group}
}
\details{
Calculation of SCORE based on the paper

,,Estimation of ten-year risk of fatal cardiovascular disease in Europe: the SCORE project''
R.M. Conroy et. al.
Eur Heart J (2003) 24 (11): 987-1003. doi: 10.1016/S0195-668X(03)00114-3
}
\author{ Przemyslaw Biecek}
\keyword{score.cardiovascular}
