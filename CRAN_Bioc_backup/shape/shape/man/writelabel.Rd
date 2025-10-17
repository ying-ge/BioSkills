\name{writelabel}
\alias{writelabel}
\title{
  adds a label next to a plot
}
\description{
  adds one-character label on left-upper margin, next to a plot
}
\usage{
writelabel(text = NULL, nr = 1, at = -0.1, line = 1, cex = 1.5, ...)
}
\arguments{
  \item{text }{text to write.
  }
  \item{nr }{integer; if text = NULL: nr is converted to uppercase letter.
  }
  \item{at }{relative distance of label position, from left margin of
    plot region.
  }
  \item{line }{line above the plot region of label position.
  }
  \item{cex }{relative size of label.
  }
  \item{... }{arguments passed to R-function \link{mtext}.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
plot(runif(2), main = "writelabel")
writelabel("A")
writelabel("B", at = 0)
writelabel("C", at = 1)}
\keyword{aplot}

