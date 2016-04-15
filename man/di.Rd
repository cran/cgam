\name{di}

\alias{di}

\title{Specify a Decreasing-Increasing Shape-Restriction in a WPS Formula}

\description{A symbolic routine to define that a surface is decreasing in one predictor and increasing in another in a formula argument to wps.}

\usage{di(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))}

\arguments{
  \item{x1}{A numeric predictor which has the same length as the response vector.}
  \item{x2}{A numeric predictor which has the same length as the response vector.}
  \item{numknots}{A vector of the number of knots used to constrain \eqn{x_1} and \eqn{x_2}. It will not be used if the user specifies the \emph{knots} argument and each predictor is within the range of its knots. The default is numknots = c(0, 0).}
  \item{knots}{A list of two vectors of knots used to constrain \eqn{x_1} and \eqn{x_2}. User-defined knots will be used if each predictor is within the range of its knots. Otherwise, \emph{numknots} and \emph{space} will be used to create knots. The default is knots = list(k1 = 0, k2 = 0).}
  \item{space}{A vector of the character specifying the method to create knots for \eqn{x_1} and \eqn{x_2}. It will not be used if the user specifies the \emph{knots} argument. If "E" is used, then equally spaced knots will be created; if "Q" is used, then a vector of equal quantiles will be created with duplicate elements removed. The number of knots is \emph{numknots} when \emph{numknots} \eqn{> 0}. Otherwise it is of the order \eqn{n^{1/6}}. The default is space = c("E", "E").}
}


\value{The vectors \eqn{x_1} and \eqn{x_2}. Each of them has six attributes, i.e., name: names of \eqn{x_1} and \eqn{x_2}; shape: "wps_di"(decreasing-increasing); numknots: the numknots argument in "di"; knots: the knots argument in "di"; space: the space argument in "di"; decreasing: two logical values indicating the monotonicity of the isotonically-constrained surface with respect to \eqn{x_1} and \eqn{x_2}, which are TRUE and FALSE.}


\details{
  "di" returns the vectors "x1" and "x2", and imposes on each vector six attributes: name, shape, numknots, knots, space and decreasing.

  The name attribute is used in the subroutine plotpersp; the numknots, knots and space attributes are the same as the numknots, knots and space arguments in "di"; the shape attribute is "wps_di"(decreasing-increasing); the decreasing values for "x1" and "x2" are TRUE and FALSE. According to the value of the vector itself and its shape, numknots, knots, space and decreasing attributes, the cone edges will be made by warped-plane spline basis functions in Meyer (2016). The cone edges are a set of basis employed in the hinge algorithm. 

  Note that "di" does not make the corresponding cone edges itself. It sets things up to a subroutine called makedelta_wps in wps.

  See references cited in this section for more details.
}


\author{
  Mary C. Meyer and Xiyue Liao
}

\references{
Meyer, M. C. (2016) Estimation and inference for isotonic regression in two dimensions, using warped-plane splines.
}

\seealso{
  \code{\link{ii}}, \code{\link{dd}}, \code{\link{wps}}
}


\examples{
  # generate data
  n <- 100
  set.seed(123)
  x1 <- sort(runif(n))
  x2 <- sort(runif(n))
  y <- 4 * (x1 - x2) + rnorm(n, sd = .5)

  # regress y on x1 and x2 under the shape-restriction: "decreasing-increasing"
  # with a penalty term = .1
  ans <- wps(y ~ di(x1, x2), pen = .1)

  # make a 3D plot of the constrained surface
  plotpersp(ans, x1, x2)
}

\keyword{shape routine}


