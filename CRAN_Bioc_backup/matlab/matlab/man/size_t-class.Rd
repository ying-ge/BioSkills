% $Id: size_t-class.Rd 23 2022-05-30 18:22:23Z proebuck $

\name{size_t-class}
\alias{size_t-class}
\docType{class}
\title{Class "size_t"}
\description{
  This class represents the dimensions of another \R object
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("size_t", \dots)}.
  Use of generator method is preferred.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{object of class \code{"integer"} containing size values}
  }
}
\section{Extends}{
Class \code{"integer"}, from data part.
Class \code{"vector"}, by class \code{"integer"}.
Class \code{"numeric"}, by class \code{"integer"}.
}
\note{
  Internal class supporting \code{\link{size}}.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{classes}

