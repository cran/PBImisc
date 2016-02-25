\name{corn}
\alias{corn}
\docType{data}
\title{A datasets relatead to gene expression in corn}
\description{
Dataset from the book ,,Modele liniowe i mieszane w R, wraz z przykladami w analizie danych''.
}
\usage{data(corn)}
\format{
data.frame with 5339 obs. and 36 variables
\describe{
   A dataset with expression of 5339 genes. Each column corresponds to a single experiment. Column name codes the setup of experiment. For example DH.C.1 is related to line DH in the condition C and it is a first technical replicate of this set of conditions.
   
   Note that a noise injection was added to this data, in order to obtain the original dataset please contact with the package maintainer.
 }
}
\details{
Dataset from the book ,,Modele liniowe i mieszane w R, wraz z przykladami w analizie danych''.

Used as an example of modeling of data from expression microarrays with the use of models with mixed effects.
}
\examples{
\dontrun{

require(lme4)

names <- colnames(corn)
X <- t(matrix(unlist(strsplit(names, ".", fixed=T)), 3, 36))
X <- data.frame(X)
colnames(X) <- c("spec", "temp", "plant")

summary(X)

y <- corn[4662,]
lmer(y~spec*temp + (1|plant:spec:temp), data=X)

}
}
\keyword{corn}

