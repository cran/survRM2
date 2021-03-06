% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.rmst2.R
\name{plot.rmst2}
\alias{plot.rmst2}
\title{plot.rmst2}
\usage{
\method{plot}{rmst2}(
  x,
  xlab = "",
  ylab = "",
  col = "red",
  col.RMST = "pink",
  col.RMTL = "orange",
  density = 80,
  angle = 85,
  ...
)
}
\arguments{
\item{x}{Results of the unadjusted analyses.}

\item{xlab}{x label.}

\item{ylab}{y label.}

\item{col}{Color for line. Default is red.}

\item{col.RMST}{Color for areas of RMST. Default is pink.}

\item{col.RMTL}{Color for areas of RMTL. Default is orange.}

\item{density}{Density of shading lines, in lines per inch. Default is 80.}

\item{angle}{Slope of shading lines, given as an angle in degrees (counter-clockwise). Default is 85.}

\item{...}{Further arguments ignored in this function.}
}
\value{
returns a plot
}
\description{
S3 method for class 'rmst2'
}
