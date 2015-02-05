\name{rmst2.sample.data}
\alias{rmst2.sample.data}
\title{Generate a sample data from the pbc data
}
\description{
The is a function to retrieve 312 randomized patients from the pbc data in survival package. 
}
\usage{
rmst2.sample.data()
}
\details{
The function creates a sample dataset to illustrate the usage of the function \code{rmst2()} in this package. 
The original pbc data in \code{survival} package consists of 418 patients data. 
This function loads the pbc data, select the 312 patients who were randomized. 
The status variable is edited, so that 1 indicates death and 0 indicates alive. 
}
\seealso{
\code{pbc} in survival package
}
\examples{
D=rmst2.sample.data()

head(D)
}
