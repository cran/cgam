\name{in.or.out}
\alias{in.or.out}

\title{To Include a Non-Parametrically Modelled Predictor in a SHAPESELECT Formula}

\description{A symbolic routine to indicate that a predictor is included as a non-parametrically modeled predictor in a formula argument to ShapeSelect.}

\usage{
in.or.out(z)
}

\arguments{
  \item{z}{A non-parametrically modelled predictor which has the same length as the response vector.}
}

\value{The vector z with three attributes, i.e., nm: the name of z; shape: 1 or 0 (in or out of the model); type: "fac" or "lin", i.e., z is modelled as a categorical predictor or a linear predictor.}

\author{Xiyue Liao}

\details{
To include a categorical predictor, in.or.out(factor(z)) is used, and to include a linear predictor z, in.or.out(z) is used. If in.or.out is not used, the user can include z in a model by adding z or factor(z) in a ShapeSelect formula.
}
\seealso{\code{\link{shapes}}, \code{\link{ShapeSelect}}}

\examples{
\dontrun{
  n <- 100
  # x is a continuous predictor 
  x <- runif(n)
  
  # generate z and to include it as a categorical predictor
  z <- rep(0:1, 50)

  # y is generated as correlated to both x and z
  # the relationship between y and x is smoothly increasing-convex
  y <- x^2 + 2 * I(z == 1) + rnorm(n, sd = 1)

  # call ShapeSelect to find the best model by the genetic algorithm
  # factor(z) may be in or out of the model  
  fit <- ShapeSelect(y ~ shapes(x) + in.or.out(factor(z)), genetic = TRUE)

  # factor(z) isn't chosen and is included in the model
  fit <- ShapeSelect(y ~ shapes(x) + factor(z), genetic = TRUE)
}
}

\keyword{variable selection routine}
