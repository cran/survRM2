\name{survRM2-package}
\alias{survRM2-package}
\docType{package}
\title{
Comparing Restricted Mean Survival Time
}

\description{
Performs two-sample comparisons using the restricted mean survival time (RMST) as a summary measure of the survival time distribution. Three kinds of between-group contrast metrics (i.e., the difference in RMST, the ratio of RMST and the ratio of the restricted mean time lost (RMTL)) are computed. The package has a function to perform an ANCOVA-type covariate adjustment as well as unadjusted analyses for those measures. 
}
\details{
\tabular{ll}{
Package: \tab survRM2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-2-5\cr
License: \tab GLP-2\cr
}
}

\author{Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui 

Maintainer: Hajime Uno <huno@jimmy.harvard.edu>
}
\references{

Uno H, Claggett B, Tian L, Inoue E, Gallo P, Miyata T, Schrag D,
Takeuchi M, Uyama Y, Zhao L, Skali H, Solomon S, Jacobus S, Hughes
M, Packer M, Wei LJ. Moving beyond the hazard ratio in quantifying the
between-group difference in survival analysis. Journal of clinical Oncology 2014, 32, 2380-2385.

Tian L, Zhao L,  Wei LJ. Predicting the restricted mean event time with
the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233.

}
\keyword{ survival }
\seealso{ survival }
\examples{ 

#--- sample data ---#
D=rmst2.sample.data()
time=D$time
status=D$status
arm=D$arm
tau=NULL
x=D[,c(4,6,7)]

#--- without covariates ----
a=rmst2(time, status, arm, tau=10)
print(a)
plot(a, xlab="Years", ylab="Probability", density=60)

#--- with covariates ----
a=rmst2(time, status, arm, tau=10, covariates=x)
print(a)

}

