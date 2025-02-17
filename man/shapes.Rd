\name{shapes}
\alias{shapes}

\title{To Include a Non-Parametrically Modelled Predictor in a SHAPESELECT Formula}

\description{A symbolic routine to indicate that a predictor is included as a non-parametrically modelled predictor in a formula argument to ShapeSelect.}

\usage{shapes(x, set = "s.9")}

\arguments{
 \item{x}{A numeric predictor which has the same length as the response vector.}
  \item{set}{A character or a numeric vector indicating all possible shapes defined for x. For example, we are not only interested in modeling the relationship between the growth of an organism (dependent variable \eqn{y}) and time (independent variable \eqn{x}), but we are also interested in the shape of the growth curve. Suppose we know \emph{a priori} that the shape could be flat, increasing, increasing concave, or increasing convex, and we further know that the curve is smooth, we can write y ~ shapes(x, set = c("flat", "s.incr", "s.incr.conc", "s.incr.conv")) in a formula to impose the four possible shape constraints on the growth curve and model it with splines. 

To be more specific, the user can choose to specify this argument as following
\describe{
\item{1. }{It could be written as "s.5", "s.9", "ord.5", "ord.9", and "tree", where "s.5" ("ord.5") means that the relationship between the response and a predictor \eqn{x} is modelled with regression splines (ordinal regression basis functions) with five possible shapes, i.e., flat, increasing, decreasing, convex, and concave; "s.9" ("ord.9") includes four more possible shapes, which are the combination of monotonicity and convexity; "tree" specifies that \eqn{x} is included as an ordinal predictor with three possibilities: no effect, tree-ordering, and unordered effect.}
\item{2. }{Or the user can choose any subset of the possible shapes, i.e., flat, increasing, decreasing, convex, concave, and combination of monotonicity and convexity. The symbols are  "flat", "incr", "decr", "conv", "conc", "incr.conv", "decr.conv", "incr.conc", and "decr.conc". To specify a spline-based regression, the user needs write something like "s.incr", "s.decr", etc.}
\item{3. }{It can also be a subset of integers between 0 and 16, where 0 is the flat shape, 1 ~ 8 indicate increasing, decreasing, convex, concave, increasing-convex, decreasing-convex, increasing-concave, and decreasing-concave, while 9 ~ 16 indicate the same shapes with a smooth assumption.}
}
}
The default is set = "s.9".
}

\value{The vector x with three attributes, i.e., nm: the name of x; shape: a numeric vector ranging from 0 to 16 to indicate possible shapes imposed on the relationship between the response and x; type: "nparam", i.e., x is non-parametrically modelled.}

\author{Xiyue Liao}

\seealso{\code{\link{in.or.out}}, \code{\link{ShapeSelect}}}

\examples{
\dontrun{
# Example 1.
  n <- 100 
   
  # generate predictors, x is non-parametrically modelled 
  # and z is parametrically modelled
  x <- runif(n)
  z <- rep(0:1, 50)
  
  # E(y) is generated as correlated to both x and z
  # the relationship between E(y) and x is smoothly increasing-convex
  y <- x^2 + 2 * I(z == 1) + rnorm(n, sd = 1)

  # call ShapeSelect to find the best model by the genetic algorithm
  fit <- ShapeSelect(y ~ shapes(x) + in.or.out(factor(z)), genetic = TRUE)

# Example 2.
  n <- 100
  z <- rep(c("A","B"), n / 2)
  x <- runif(n)

  # y0 is generated as correlated to z with a tree-ordering in it
  # y0 is smoothly increasing-convex in x
  y0 <- x^2 + I(z == "B") * 1.5
  y <- y0 + rnorm(n, 1)

  fit <- ShapeSelect(y ~ s.incr(x) + shapes(z, set = "tree"), genetic = FALSE)
  
  # check the best fit in terms of z
  fit$top
}
}

\keyword{variable selection routine}
