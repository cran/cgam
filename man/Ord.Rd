\name{Ord}
\alias{Ord}
\title{Specify an Ordered Categorical Family in a CGAM Formula}

\description{
This is a subroutine to specify an ordered catergorical family in a cgam formula. It set things up to a routine called cgam.polr. This is learned from the polr routine in the MASS package, which fits a logistic or probit regression model to an ordered categorical response. Currently only the logistic regression model is allowed. 
}

\usage{Ord(link = "identity")}

\arguments{
  \item{link}{The link function. Users don't need specify this term.}
}

\details{
See the polr section in the official manual of the MASS package (\url{https://cran.r-project.org/package=MASS}) for details. 
}

\value{
\item{muhat}{The estimated expected value of a latent variable.}
\item{zeta}{Estimated cut-points defining the intervals of a latent 
variable such that the latent variable is between two adjacent cut-points is equivalent to that the ordered categorical response is in a category.}
}

\author{
Xiyue Liao
}

\references{
	Agresti, A. (2002) \emph{Categorical Data}. Second edition. Wiley.
}


\seealso{
  \code{\link{mental}}
}

\examples{
\dontrun{
	# Example 1. 
	# generate the predictor and the latenet variable
	n <- 500
	set.seed(123)
	x <- runif(n, 0, 1)
	yst <- 5*x^2 + rlogis(n)
	
	# generate observed ordered response, which has levels 1, 2, 3. 
	cts <- quantile(yst, probs = seq(0, 1, length = 4)) 
	yord <- cut(yst, breaks = cts, include.lowest = TRUE, labels = c(1:3), Ord = TRUE)
	y <- as.numeric(levels(yord))[yord] 

	# regress y on x under the shape-restriction: the latent variable is "increasing-convex"
	# w.r.t x
	ans <- cgam(y ~ s.incr.conv(x), family = Ord)
	
	# check the estimated cut-points
	ans$zeta
	
	# check the estimated expected value of the latent variable
	head(ans$muhat)
	
	# check the estimated probabilities P(y = k), k = 1, 2, 3
	head(fitted(ans))
	
	# check the estimated latent variable
	plot(x, yst, cex = 1, type = "n", ylab = "latent variable")
	cols <- topo.colors(3)
	for (i in 1:3) {
		points(x[y == i], yst[y == i], col = cols[i], pch = i, cex = 0.7)
	}
	for (i in 1:2) {
		abline(h = (ans$zeta)[i], lty = 4, lwd = 1)
	}
	lines(sort(x), (5*x^2)[order(x)], lwd = 2)
	lines(sort(x), (ans$muhat)[order(x)], col = 2, lty = 2, lwd = 2)
	legend("topleft", bty = "n", col = c(1, 2), lty = c(1, 2), 
	c("true latent variable", "increasing-convex fit"), lwd = c(1, 1))
}	
\dontrun{
	# Example 2. mental impairment data set 
	# mental impairment is an ordinal response with 4 categories recorded as 1, 2, 3, and 4
	# two covariates are life event index and socio-economic status (high = 1, low = 0)
	data(mental)
	table(mental$mental)
	
	# model the relationship between the latent variable and life event index as increasing
	# socio-economic status is included as a binary covariate 
	fit.incr <- cgam(mental ~ incr(life) + ses, data = mental, family = Ord)

	# check the estimated probabilities P(mental = k), k = 1, 2, 3, 4
	probs.incr <- fitted(fit.incr)
	head(probs.incr)
}

}

\keyword{ordered categorical family}