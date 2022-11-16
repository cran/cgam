\name{predict.cgam}
\alias{predict.cgam}

\title{Predict Method for CGAM Fits}

\description{Predicted values based on a cgam object}

\usage{
\method{predict}{cgam}(object, newData, interval = c("none", "confidence", "prediction"), 
type = c("response", "link"), level = 0.95, n.mix = 500,...)
}

\arguments{
  \item{object}{A cgam object.}
  \item{newData}{A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{interval}{Type of interval calculation. A prediction interval is only implemented for Gaussian response for now.}
  \item{type}{If the response is Gaussian, type = "response" gives the predicted mean; if the response is binomial, type = "response" gives the predicted probabilities, and type = "link" gives the predicted systematic component.}
  \item{level}{Tolerance/confidence level.}
  \item{n.mix}{Number of simulations to get the mixture distribution. The default is n.mix = 500.}
   \item{\dots}{Further arguments passed to the routine.}
}

\details{
	Constrained spline estimators can be characterized by projections onto a polyhedral convex cone. Point-wise confidence intervals for constrained splines are constructed by estimating the probabilities that the projection lands on each of the faces of the cone, and using a mixture of covariance matrices to estimate the standard error of the function estimator at any design point. 
	
	Note that currently predict.cgam only works when all components in a cgam formula are additive.

	See references cited in this section for more details.	
}

\value{
\item{fit}{A vector of predictions.}
\item{lower}{A vector of lower bound if interval is set to be "confidence".}
\item{upper}{A vector of upper bound if interval is set to be "confidence".}
}

\references{
	Meyer, M. C. (2017) Constrained partial linear regression splines. \emph{Statistical Sinica in press}.

	Meyer, M. C. (2017) Confidence intervals for regression functions using constrained
splines with application to estimation of tree height

	Meyer, M. C. (2012) Constrained penalized splines. \emph{Canadian Journal of Statistics \bold{40(1)}}, 190--206.
	  	  
	Meyer, M. C. (2008) Inference using shape-restricted regression splines. \emph{Annals of Applied Statistics \bold{2(3)}}, 1013--1033.
}

\author{
 Mary C. Meyer and Xiyue Liao
}

\examples{
# Example 1.
	# generate data
	n <- 100
	set.seed(123)
	x <- runif(n)
	y <- 4*x^3 + rnorm(n)

	# regress y on x under the shape-restriction: "increasing-convex"
	fit <- cgam(y ~ s.incr.conv(x))
	
	# make a data frame 
	x0 <- seq(min(x), max(x), by = 0.05)
	new.Data <- data.frame(x = x0)

	# predict values in new.Data based on the cgam fit without a confidence interval
	pfit <- predict(fit, new.Data)
	
	# or 
	pfit <- predict(fit, new.Data, interval = "none")
	
	# make a plot to check the prediction
	plot(x, y, main = "Predict Method for CGAM")
	lines(sort(x), (fitted(fit)[order(x)]))
	points(x0, pfit$fit, col = 2, pch = 20)

	# predict values in new.Data based on the cgam fit with a 95 percent confidence interval
	pfit <- predict(fit, new.Data, interval = "confidence", level = 0.95)

	# make a plot to check the prediction
	plot(x, y, main = "Pointwise Confidence Bands (Gaussian Response)")
	lines(sort(x), (fitted(fit)[order(x)]))
	lines(sort(x0), (pfit$lower)[order(x0)], col = 2, lty = 2)
	lines(sort(x0), (pfit$upper)[order(x0)], col = 2, lty = 2)
	points(x0, pfit$fit, col = 2, pch = 20)

# Example 2. binomial response
	n <- 200
	x <- seq(0, 1, length = n)

	eta <- 4*x - 2
	mu <- exp(eta)/(1+exp(eta))
	set.seed(123)
	y <- 1:n*0
	y[runif(n)<mu] = 1
        
	fit <- cgam(y ~ s.incr.conv(x), family = binomial)
	muhat <- fitted(fit)
	
    # predict values in new.Data based on the cgam fit with a 95 percent confidence interval
    xinterp <- seq(min(x), max(x), by = 0.05)
    new.Data <- data.frame(x = xinterp)
	pfit <- predict(fit, new.Data, interval = "confidence", level = 0.95)
	pmu <- pfit$fit
	lwr <- pfit$lower
	upp <- pfit$upper
	
	# make a plot to check the prediction	
	plot(x, y, type = "n", ylim = c(0, 1), 
	main = "Pointwise Confidence Bands (Binomial Response)")
	rug(x[y == 0])
    rug(x[y == 1], side = 3)   
	lines(x, mu)
	lines(x, muhat, col = 5, lty = 2)
       
	points(xinterp, pmu, pch = 20)
	lines(xinterp, upp, col = 5)
	points(xinterp, upp, pch = 20)
	lines(xinterp, lwr, col = 5)
	points(xinterp, lwr, pch = 20)
}

\keyword{prediction routine}