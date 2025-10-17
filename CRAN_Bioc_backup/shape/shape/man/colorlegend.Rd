\name{colorlegend}
\alias{colorlegend}
\title{
  adds a color legend to a plot.
}
\description{
  Adds a color legend to a plot.
}
\usage{
colorlegend(col = femmecol(100), zlim, zlevels = 5, dz = NULL,
            zval = NULL, log = FALSE, posx = c(0.9, 0.93), 
            posy = c(0.05, 0.9), main = NULL, main.cex = 1.0, 
            main.col = "black", lab.col = "black", 
            digit = 0, left = FALSE, ...)
}

\arguments{
  \item{col }{color palette to be used; also allowed are two extremes
    or one value.
  }
  \item{zlim }{two-valued vector, the minimum and maximum z values.
  }
  \item{zlevels }{number of z-levels, one value, ignored if \code{dz} or
    \code{zval} not equal to \code{NULL}.
  }
  \item{dz }{increment in legend values, one value; ignored if \code{zval}
    not equal to \code{NULL}.
  }
  \item{zval }{a vector of z-values to label legend.
  }
  \item{log }{logical indicating whether to log transform or not.
  }
  \item{posx }{relative position of left and right edge of color bar on
    first axis, [0,1].
  }
  \item{posy }{relative position on lower and upper edge of colar bar on
    second axis, [0,1].
  }
  \item{main }{main title, written above the color bar.
  }
  \item{main.cex }{relative size of main title.
  }
  \item{main.col }{color of main title.
  }
  \item{lab.col }{color of labels.
  }
  \item{digit }{number of significant digits in labels.
  }
  \item{left }{logical indicating whether to put the labels on the right
    (\code{TRUE}) or on the left (\code{FALSE}).
  }
  \item{... }{arguments passed to R-function \link{text} when writing labels.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
emptyplot(main = "colorlegend")
colorlegend(zlim = c(0, 10))
colorlegend(posx = c(0.8, 0.83), col = greycol(100), 
            zlim = c(0, 1), digit = 1)
colorlegend(posx = c(0.7, 0.73), left = TRUE, col = rainbow(100), 
            zlim = c(0, 10), digit = 1, dz = 2.5)
colorlegend(posx = c(0.5, 0.53), 
            col = intpalette(c("red", "yellow", "black"), 100), 
            zlim = c(0, 20), zval = c(1, 3, 7, 15))
colorlegend(posy = c(0.0, 0.15), posx = c(0.2, 0.3), 
            col = rainbow(100), zlim = c(0, 1), 
            zlevels = NULL, main = "rainbow")
colorlegend(posy = c(0.25, 0.4), posx = c(0.2, 0.3), 
            zlim = c(0, 1), zlevels = NULL, main = "femmecol")
colorlegend(posy = c(0.5, 0.65), posx = c(0.2, 0.3), 
            col = terrain.colors(100), zlim = c(0, 1), 
            zlevels = NULL, main = "terrain.colors")
colorlegend(posy = c(0.75, 0.9), posx = c(0.2, 0.3), 
            col = heat.colors(100), zlim = c(0, 1), 
            zlevels = NULL, main = "heat.colors")
}

\keyword{aplot}

