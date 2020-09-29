\name{s.decr.decr}

\alias{s.decr.decr}

\title{Specify a Doubly-Decreasing Shape-Restriction in a CGAM Formula}

\description{A symbolic routine to define that a surface is decreasing in two predictors in a formula argument to cgam.}

\usage{s.decr.decr(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))}

\arguments{
  \item{x1}{A numeric predictor which has the same length as the response vector.}
  \item{x2}{A numeric predictor which has the same length as the response vector.}
  \item{numknots}{A vector of the number of knots used to constrain \eqn{x_1} and \eqn{x_2}. It will not be used if the user specifies the \emph{knots} argument and each predictor is within the range of its knots. The default is numknots = c(0, 0).}
  \item{knots}{A list of two vectors of knots used to constrain \eqn{x_1} and \eqn{x_2}. User-defined knots will be used if each predictor is within the range of its knots. Otherwise, \emph{numknots} and \emph{space} will be used to create knots. The default is knots = list(k1 = 0, k2 = 0).}
  \item{space}{A vector of the character specifying the method to create knots for \eqn{x_1} and \eqn{x_2}. It will not be used if the user specifies the \emph{knots} argument. If "E" is used, then equally spaced knots will be created; if "Q" is used, then a vector of equal quantiles will be created with duplicate elements removed. The number of knots is \emph{numknots} when \emph{numknots} is a positive integer \eqn{> 4}. Otherwise it is of the order \eqn{n^{1/6}}. The default is space = c("E", "E").}
}


\value{The vectors \eqn{x_1} and \eqn{x_2}. Each of them has six attributes, i.e., name: names of \eqn{x_1} and \eqn{x_2}; shape: "wps_dd"(doubly-decreasing); numknots: the numknots argument in "s.decr.decr"; knots: the knots argument in "s.decr.decr"; space: the space argument in "s.decr.decr"; decreasing: two logical values indicating the monotonicity of the isotonically-constrained surface with respect to \eqn{x_1} and \eqn{x_2}, which are both TRUE.}


\details{
  "s.decr.decr" returns the vectors "x1" and "x2", and imposes on each vector six attributes: name, shape, numknots, knots, space and decreasing.

  The name attribute is used in the subroutine plotpersp; the numknots, knots and space attributes are the same as the numknots, knots and space arguments in "s.decr.decr"; the shape attribute is "wps_dd"(doubly-decreasing); the decreasing values for both vectors are TRUE. According to the value of the vector itself and its shape, numknots, knots, space and decreasing attributes, the cone edges will be made by warped-plane spline basis functions in Meyer (2016). The cone edges are a set of basis employed in the hinge algorithm. 

  Note that "s.decr.decr" does not make the corresponding cone edges itself. It sets things up to a subroutine called makedelta_wps.

  See references cited in this section for more details.
}


\author{
  Mary C. Meyer and Xiyue Liao
}

\references{
Meyer, M. C. (2017) Estimation and inference for regression surfaces using shape-constrained splines.}

\seealso{
  \code{\link{s.incr.incr}}, \code{\link{s.decr.incr}}, \code{\link{s.incr.decr}}, \code{\link{cgam}}
}


\examples{
\dontrun{
  # generate data
  n <- 100
  set.seed(123)
  x1 <- runif(n)
  x2 <- runif(n)
  y <- -4 * (x1 + x2 - x1 * x2) + rnorm(n, sd = .2)

  # regress y on x1 and x2 under the shape-restriction: "doubly-decreasing" 
  # using the penalized estimator
  ans <- cgam(y ~ s.decr.decr(x1, x2), pnt = TRUE)

  # make a 3D plot of the constrained surface
  plotpersp(ans)	
}
}

\keyword{shape routine}


