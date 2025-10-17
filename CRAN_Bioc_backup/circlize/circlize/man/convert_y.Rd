\name{convert_y}
\alias{convert_y}
\title{
Convert unit on y direction in data coordinate
}
\description{
Convert unit on y direction in data coordinate
}
\usage{
convert_y(
    x,
    unit = c("mm", "cm", "inches"),
    sector.index = get.current.sector.index(),
    track.index = get.current.track.index())
}
\arguments{

  \item{x}{a numeric vector}
  \item{unit}{supported units, only "mm", "cm", "inches"}
  \item{sector.index}{index for the sector where the conversion is applied}
  \item{track.index}{index for the track where the conversion is applied}

}
\value{
A vector of numeric values which are measured in the specified data coordinate
}
\seealso{
For pre-defined units, users can use \code{\link{cm_y}}, \code{\link{mm_y}} and \code{\link{inches_y}}.

\code{\link{convert_x}} converts on x direction.

\url{https://jokergoo.github.io/circlize_book/book/circular-layout.html#convert-functions}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see example on `convert_x` page 
NULL
}
