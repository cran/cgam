\name{s.decr}

\alias{s.decr}

\title{Specify a Smooth and Decreasing Shape-Restriction in a CGAM Formula}

\description{A symbolic routine to define that the systematic component \eqn{\eta} is smooth and decreasing in a predictor in a formula argument to cgam. This is the smooth version.}

\usage{s.decr(x, numknots = 0, knots = 0, var.knots = 0, space = "E", db.exp = FALSE)}

\arguments{
  \item{x}{A numeric predictor which has the same length as the response vector.}
 \item{numknots}{The number of knots used to constrain \eqn{x}. It will not be used if the user specifies the \emph{knots} argument. The default is numknots = \eqn{0}.}
  \item{knots}{The knots used to constrain \eqn{x}. User-defined knots will be used when given. Otherwise, \emph{numknots} and \emph{space} will be used to create knots. The default is knots = \eqn{0}.}
    \item{var.knots}{The knots used in variance function estimation. User-defined knots will be used when given. The default is var.knots = \eqn{0}.} 
 \item{space}{A character specifying the method to create knots. It will not be used if the user specifies the \emph{knots} argument. If space == "E", then equally spaced knots will be created; if space == "Q", then a vector of equal \eqn{x} quantiles will be created based on \eqn{x} with duplicate elements removed. The number of knots is \emph{numknots} when \emph{numknots} \eqn{> 0}. Otherwise it is of the order \eqn{n^{1/7}}. The default is space = "E".}
\item{db.exp}{The parameter will be used in variance function estimation. If db.exp = TRUE, then the errors are assumed to follow a normal distribution; otherwise, the errors are assumed to follow a double-exponential distribution. The default is db.exp = FALSE.}
}

\value{The vector x with five attributes, i.e, name: the name of x; shape: 10("smooth and decreasing"); numknots: the numknots argument in "s.decr"; knots: the knots argument in "s.decr"; space: the space argument in "s.decr".}

\details{
  "s.decr" returns the vector "x" and imposes on it five attributes: name, shape, numknots, knots, space, var.knots and db.exp. 

  The name attribute is used in the subroutine plotpersp; the numknots, knots and space attributes are the same as the numknots, knots and space arguments in "s.decr"; the shape attribute is 10("smooth and decreasing"). According to the value of the vector itself and its shape, numknots, knots and space attributes, the cone edges will be made by I-spline basis functions in Meyer (2008). The cone edges are a set of basis employed in the hinge algorithm. 

  Note that "s.decr" does not make the corresponding cone edges itself. It sets things up to a subroutine called makedelta in cgam.

  var.knots and db.exp will be used for monotonic variance function estimation.
  
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

  # regress y on x under the shape-restriction: "smooth and decreasing"
  ans <- cgam(y ~ s.decr(x))
  knots <- ans$knots[[1]]

  # make a plot
  par(mar = c(4, 4, 1, 1))
  plot(x, y, cex = .7, xlab = "x", ylab = "y")
  lines(x, ans$muhat, col = 2)
  legend("topleft", bty = "n", "smooth and decreasing fit", col = 2, lty = 1)
  legend(-.3, 8, bty = "o", "knots", pch = "X")
  points(knots, 1:length(knots)*0+min(y), pch = "X")
}

\seealso{
  \code{\link{decr}}
}

\keyword{shape routine}
