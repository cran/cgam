\name{s.decr.conv}
\alias{s.decr.conv}


\title{Specify a Smooth, Decreasing and Convex Shape-Restriction in a CGAM Formula}

\description{A symbolic routine to define that the systematic component \eqn{\eta} is smooth, decreasing and convex in a predictor in a formula argument to cgam. This is the smooth version.}

\usage{s.decr.conv(x, numknots = 0, knots = 0, space = "Q")}


\arguments{
  \item{x}{A numeric predictor which has the same length as the response vector.}
 \item{numknots}{The number of knots used to constrain \eqn{x}. It will not be used if the user specifies the \emph{knots} argument. The default is numknots = \eqn{0}.}
  \item{knots}{The knots used to constrain \eqn{x}. User-defined knots will be used when given. Otherwise, \emph{numknots} and \emph{space} will be used to create knots. The default is knots = \eqn{0}.}
 \item{space}{A character specifying the method to create knots. It will not be used if the user specifies the \emph{knots} argument. If space == "E", then equally spaced knots will be created; if space == "Q", then a vector of equal \eqn{x} quantiles will be created based on \eqn{x} with duplicate elements removed. The number of knots is \emph{numknots} when \emph{numknots} \eqn{> 0}. Otherwise it is of the order \eqn{n^{1/7}}. The default is space = "Q".}
}

\value{The vector x with five attributes, i.e., name: the name of x; shape: 15("smooth, decreasing and convex"); numknots: the numknots argument in "s.decr.conv"; knots: the knots argument in "s.decr.conv"; space: the space argument in "s.decr.conv".}

\details{
 "s.decr.conv" returns the vector "x" and imposes on it five attributes: name, shape, numknots, knots and space.

  The name attribute is used in the subroutine plotpersp; the numknots, knots and space attributes are the same as the numknots, knots and space arguments in "s.decr.conv"; the shape attribute is 15("smooth, decreasing and convex"). According to the value of the vector itself and its shape, numknots, knots and space attributes, the cone edges will be made by C-spline basis functions in Meyer (2008). The cone edges are a set of basis employed in the hinge algorithm.

  Note that "s.decr.conv" does not make the corresponding cone edges itself. It sets things up to a subroutine called makedelta in cgam.

  See references cited in this section for more details.
}


\author{
  Mary C. Meyer and Xiyue Liao
}

\references{
   Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications in Statistics \bold{42(5)}}, 1126--1139.

  Meyer, M. C. (2008) Inference using shape-restricted regression splines. \emph{Annals of Applied Statistics \bold{2(3)}}, 1013--1033.
}

\examples{
  data(cubic)

  # extract x
  x <- - cubic$x

  # extract y
  y <- cubic$y

  # regress y on x under the shape-restriction: "smooth, decreasing and convex"
  ans <- cgam(y ~ s.decr.conv(x))
  knots <- ans$knots[[1]]

  # make a plot
  par(mar = c(4, 4, 1, 1))
  plot(x, y, cex = .7, xlab = "x", ylab = "y")
  lines(x, ans$muhat, col = 2)
  legend("topleft", bty = "n", "smooth, decreasing and convex fit", col = 2, lty = 1)
  legend(-.3, 9.2, bty = "o", "knots", pch = "X")
  points(knots, 1:length(knots)*0+min(y), pch = "X")
}

\seealso{
  \code{\link{decr.conv}}
}

\keyword{shape routine}
