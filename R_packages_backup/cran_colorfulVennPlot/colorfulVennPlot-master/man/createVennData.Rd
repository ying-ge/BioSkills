\name{createVennData}
\alias{createVennData}
\title{Generates data for plotVenn3d()}
\description{Given data with 3 columns, it will generate a vector of 7 numerical values, which describe grouping of the data across 3 specified columns. The output should be used for plotVenn3d() function.
}
\usage{
createVennData(x, Cols = NULL, Splits = c(0, 0, 0), Labels = NULL, 
  type = c("count", "percent"), ToSkip = "000")
}
\arguments{
  \item{x}{ a table of numerical values of at least 3 columns. }
  \item{Cols}{ a numeric or character vector of length 3, specifying column names or numbers. If not specified, it defaults to 3 first columns of input data 'x'. }
  \item{Splits}{ a numeric vector of length 3, specifying splitting value for each column. See below for more information. }
  \item{Labels}{ a character vector of length 3, clarifying names of each group/column. If not specified, it defaults to value of variable 'Cols'. }
  \item{type}{ specifies the type of output data. 
  'count' counts a number of each grouping
  'percent' computes the percents of each grouping}
  \item{ToSkip}{ defaults to '000'. Do not change it, unless you know exactly what you are doing. }
}
\details{
NOTICE: This only works for 3-dimensional data.
Data 'x' should be a table of numerical values of at least 3 columns.

Each column's data is split by corresponding 'Splits' value into 'above' and 'below'. Then, it is gathered together in the form:

  1 - how many rows of data 'x' had all 3 columns 'above' ("111")

  2 - how many rows of data 'x' had 1st column 'below' Splits[1] and columns 2 & 3 'above' Splits[2:3] ("011")

  3 - how many rows of data 'x' had 1st column 'above' Splits[1], 2nd 'below' Splits[2] and 3rd 'above' Splits[3] ("101") 

  4 - how many rows of data 'x' had 1st column 'below' Splits[1], 2nd 'below' Splits[2] and 3rd 'above' Splits[3] ("001")

  5 - how many rows of data 'x' had 1st column 'above' Splits[1], 2nd 'above' Splits[2] and 3rd 'below' Splits[3] ("110")

  6 - how many rows of data 'x' had 1st column 'below' Splits[1], 2nd 'above' Splits[2] and 3rd 'below' Splits[3] ("010")

  7 - how many rows of data 'x' had 1st column 'above' Splits[1], 2nd 'below' Splits[2] and 3rd 'below' Splits[3] ("100")

  8 - how many rows of data 'x' had all 3 columns 'below' ("000")

Once this vector is generated, it transforms it into percents, if specified in variable 'type'. At this point, the sum will add up to 100 percent.

Then, one value has to be dropped - by definition, a 3-dimensional Venn diagram plots 7 out of 8 values, the space around the plot being the last value. Generally, value '000' must be dropped - it is set as default, and we do not recommend changing it.

The output will be then a vector of length 7 rather than 8, missing the 'ToSkip' part.
}
\value{
  a list containing:
  
  \item{x}{ a numeric vector of length 7 describing the counts/percents of split data. See 'Details' for more information. }
  \item{labels}{ a character vector of 3, specifying the names of the groups. }
}
\author{
  Aliona Manvae
}
\examples{
  data(survey, package = "MASS")
  # For more info: ?MASS::survey
  
  # Grouping students by Pulse, Height and Age
  plot.new()
  vennData <- createVennData(survey,
    Cols = c("Pulse","Height","Age"),
    Splits = c(75, 170, 20),   ## Split Pulse at 75, split Height at 170, split Age at 20
    Labels = NULL,   ## Keep Labels NULL, as we wish to stick with Pulse, Height and Age titles.
    type = c("percent"))   ## generate data as percents

  # Plot 
  plot.new()
  plotVenn3d(vennData$x, vennData$labels)
  # OR:
  plot.new()
  plotVenn(vennData$x, vennData$labels)
  
}
\keyword{ manip }
