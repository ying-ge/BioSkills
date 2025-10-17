\name{convert_length}
\alias{convert_length}
\title{
Convert units
}
\description{
Convert units
}
\usage{
convert_length(x, unit = c("mm", "cm", "inches"))
}
\arguments{

  \item{x}{a numeric vector}
  \item{unit}{supported units, only "mm", "cm", "inches".}

}
\details{
This function coverts mm/cm/inches units to units measured in the canvas coordinate,
e.g. how much is it in the canvas coordinate for 1 mm/cm/inches.

Since in the circular plot, the aspect ratio is always 1, it does not matter this conversion
is applied on x direction or y direction.

This function is mainly used in the radical direction.
}
\seealso{
\code{\link{convert_x}} and \code{\link{convert_y}} convert absolute units into a data coordinate in a specified cell.

\url{https://jokergoo.github.io/circlize_book/book/circular-layout.html#convert-functions}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
sectors = letters[1:10]
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
circos.initialize(sectors, xlim = cbind(rep(0, 10), runif(10, 0.5, 1.5)))
circos.track(ylim = c(0, 1), track.height = mm_h(5))
circos.par(track.margin = c(0, mm_h(2)))
circos.track(ylim = c(0, 1), track.height = cm_h(1))
circos.par(track.margin = c(0, mm_h(5)))
circos.track(ylim = c(0, 1), track.height = inch_h(1))
circos.clear()
}
