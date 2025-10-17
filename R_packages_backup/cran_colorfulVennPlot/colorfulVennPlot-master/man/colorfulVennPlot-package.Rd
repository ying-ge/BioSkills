\name{colorfulVennPlot-package}
\alias{colorfulVennPlot-package}
\alias{colorfulVennPlot}
\docType{package}
\title{ Venn diagram with user-specified colors for 2-dimensional, 3-dimensional and 4-dimensional data }
\description{ Plot Venn diagram of 2-dimensional, 3-dimensional and 4-dimensional data with user-specified custom colors. }
\details{
\tabular{ll}{
Package: \tab plotVenn\cr
Type: \tab Package\cr
Version: \tab 2.4\cr
Date: \tab 2013-11-11\cr
License: \tab GPL-2\cr
}
The package is designed to be able to plot Venn diagrams for 2-d, 3-d and 4-d data. The uniqueness of it is in being able to specify custom colors for each part of the plot - including the crossovers. 

If data needs to be shaped into appropriate format for it, createVennData() can be used. Given 3 columns of data and values by which to split the data, createVennData() splits it into 7 groups and counts values in each group or percent allocations. Notice: currently, it only can be used to generate 3-dimensional data.

}
\author{
Elliot Noma, Aliona Manvae

Maintainer: Elliot Noma <noma@garrettassetmanagement.com>
}
\keyword{ package }
\examples{
  data(survey, package = "MASS")
  # For more info: ?MASS::survey
  
  # Grouping students by Pulse, Height and Age
  vennData <- createVennData(survey,
    Cols = c("Pulse","Height","Age"),
    Splits = c(75, 170, 20),   ## Split Pulse at 75, split Height at 170, split Age at 20
    Labels = NULL,   ## Keep Labels NULL, as we wish to stick with Pulse, Height and Age titles.
    type = c("percent"))   ## generate data as percents

  # Plot Venn diagram
  # Notice: 
  # plotVenn is a generic function, which can be used to call plotVenn2d, plotVenn3d or plotVenn4d
  plot.new()
  plotVenn(vennData$x, vennData$labels)
  
  ##############################################
  # Example of 4-dimensional data plotting, with default data
  plot.new()
  plotVenn4d(Title = "Example of 4-d Venn diagram")

}
