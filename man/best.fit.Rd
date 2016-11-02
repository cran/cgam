\name{best.fit}
\alias{best.fit}
\title{Extract the Best Fit Returned by the ShapeSelect Routine}

\description{The is a subroutine which only works for the ShapeSelect routine. It returns an object of the cgam class given the variables and their shapes chosen by the ShapeSelect routine.
}

\usage{best.fit(x)}

\arguments{
  \item{x}{x is an object of the ShapeSelect class.}
}

\value{
  \item{object}{The best fit returned by the ShapeSelect routine, which is an object of the cgam class.}
}


\author{
Xiyue Liao
}

\examples{
\dontrun{
  library(MASS)
  data(Rubber)

  # do a variable and shape selection with four possible shapes
  # increasing, decreasing, convex and concave 
  ans <- ShapeSelect(loss ~ shapes(hard, set = c("incr", "decr", "conv", "conc"))
  + shapes(tens, set = c("incr", "decr", "conv", "conc")), data = Rubber, genetic = TRUE)

  # check the best fit, which is an object of the cgam class
  bf <- best.fit(ans)
  class(bf)
  plotpersp(bf)
}
}

\seealso{\code{\link{cgam}}, \code{\link{ShapeSelect}}}

\keyword{best fit of the ShapeSelect routine}
