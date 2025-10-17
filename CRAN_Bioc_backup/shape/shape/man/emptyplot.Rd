\name{emptyplot}
\alias{emptyplot}
\title{open a plot without axes, labels,...}
\description{
  Creates a plotting region, bounded by xlim and ylim;
  without axes, labels, titles, useful for plotting shapes.
}
\usage{
emptyplot(xlim = c(0, 1), ylim = xlim, asp = 1, frame.plot = FALSE,
          col = NULL, ...)
}

\arguments{
  \item{xlim }{the x limits (min,max) of the plot.
  }
  \item{ylim }{the y limits (min,max) of the plot.
  }
  \item{asp }{the y/x aspect ratio.
  }
  \item{frame.plot }{to toggle off drawing of a bounding box.
  }
  \item{col }{the background color.
  }
  \item{... }{arguments passed to R-function \link{plot}.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\seealso{
  \code{\link{plot}}, \code{\link{plot.default}}
}

\keyword{hplot}

