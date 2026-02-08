\name{cgamm}

\alias{cgamm}

\title{Constrained Generalized Additive Mixed-Effects Model Fitting}

\description{
	This routine is an addition to the main routine cgam in this package. A random-intercept component is included in a cgam model.
}

\usage{
cgamm(formula, nsim = 100, family = gaussian(), cpar = 1.2, data = NULL, weights = NULL,
sc_x = FALSE, sc_y = FALSE, bisect = TRUE, reml = TRUE, 
multicore = getOption("cgam.parallel"), CIC = FALSE, nAGQ = 1L)
}

\arguments{
  \item{formula}{A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor + (1|id)", where id is the label for a group effect. For now, only gaussian responses are considered and this routine only includes a random-intercept effect. See \code{\link{cgam}} for more details.}
  \item{nsim}{The number of simulations used to get the cic parameter. The default is nsim = 0.}
  \item{family}{A parameter indicating the error distribution and link function to be used in the model. For now, the options are family = gaussian() and family = binomial().}
  \item{cpar}{A multiplier to estimate the model variance, which is defined as \eqn{\sigma^2 = SSR / (n - cpar * edf)}. SSR is the sum of squared residuals for the full model and edf is the effective degrees of freedom. The default is cpar = 1.2. The user-defined value must be between 1 and 2. See Meyer, M. C. and M. Woodroofe (2000) for more details.}
  \item{data}{An optional data frame, list or environment containing the variables in the model. The default is data = NULL.}
  \item{weights}{An optional non-negative vector of "replicate weights" which has the same length as the response vector. If weights are not given, all weights are taken to equal 1. The default is weights = NULL.}
  \item{sc_x}{Logical flag indicating if or not continuous predictors are normalized. The default is sc_x = FALSE.}
  \item{sc_y}{Logical flag indicating if or not the response variable is normalized. The default is sc_y = FALSE.}
  \item{bisect}{If bisect = TRUE, a 95 percent confidence interval will be found for the variance ratio parameter by a bisection method.}
  \item{reml}{If reml = TRUE, restricted maximum likelihood (REML) method will be used to find estimates instead of maximum likelihood estimation (MLE).}
  \item{multicore}{
A parameter retrieving the current global option for "csurvey.multicore" and assigns its value to the multicore variable, allowing the function to respect a user-defined setting for parallel processing behavior across the package.
}
\item{CIC}{
Logical flag indicating if or not cic is computed. 
}
  \item{nAGQ}{Integer scalar - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. Defaults to 1, corresponding to the Laplace approximation. Values greater than 1 produce greater accuracy in the evaluation of the log-likelihood at the expense of speed.}
}

\value{
  \item{muhat}{The fitted fixed-effect term.}
  \item{ahat}{A vector of estimated random-effect terms.}
  \item{sig2hat}{Estimate of the variance (\eqn{\sigma^2}) of between-cluster error terms.}
  \item{siga2hat}{Estimate of the variance (\eqn{\sigma_a^2}) of within-cluster error terms.}
  \item{thhat}{Estimate of the ratio (\eqn{\theta}) of two variances.}
  \item{pv.siga2}{\eqn{p}-value of the test \eqn{H_0: \sigma_a^2=0}}
  \item{ci.siga2}{95 percent confidence interval for the variance of within-cluster error terms.}
  \item{ci.th}{95 percent confidence interval for ratio of two variances.}
  \item{ci.rho}{95 percent confidence interval for intraclass correlation coefficient.}
  \item{ci.sig2}{95 percent confidence interval for the variance of between-cluster error terms.}
  \item{call}{The matched call.}
}

\examples{
# Example 1 (family = gaussian).

# simulate a balanced data set with 30 clusters
# each cluster has 30 data points
	n <- 30
	m <- 30

# the standard deviation of between cluster error terms is 1
# the standard deviation of within cluster error terms is 2
	sige <- 1
	siga <- 2

# generate a continuous predictor
	x <- 1:(m*n)
	for(i in 1:m) {
		x[(n*(i-1)+1):(n*i)] <- round(runif(n), 3)
	}
# generate a group factor
	group <- trunc(0:((m*n)-1)/n)+1

# generate the fixed-effect term
	mu <- 10*exp(10*x-5)/(1+exp(10*x-5))

# generate the random-intercept term asscosiated with each group
	avals <- rnorm(m, 0, siga)

# generate the response
	y <- 1:(m*n)
	for(i in 1:m){
		y[group == i] <- mu[group == i] + avals[i] + rnorm(n, 0, sige)
	}

# use REML method to fit the model
	ans <- cgamm(y ~ s.incr(x) + (1|group), reml=TRUE)
	summary(ans)
	anova(ans)

	muhat <- ans$muhat
	plot(x, y, col = group, cex = .6)
	lines(sort(x), mu[order(x)], lwd = 2)
	lines(sort(x), muhat[order(x)], col = 2, lty = 2, lwd = 2)
	legend("topleft", bty = "n", c("true fixed-effect term", "cgamm fit"),
	col = c(1, 2), lty = c(1, 2), lwd = c(2, 2))

# Example 2 (family = binomial).
# simulate a balanced data set with 20 clusters
# each cluster has 20 data points

  n <- 20
  m <- 20#
  N <- n*m

  # siga is the sd for the random intercept
  siga <- 1

# generate a group factor
  group <- trunc(0:((m*n)-1)/n)+1
  group <- factor(group)

# generate the random-intercept term asscosiated with each group
  avals <- rnorm(m,0,siga)

# generate the fixed-effect mean term: mu, systematic term: eta and the response: y
  x <- runif(m*n)
  mu <- 1:(m*n)
  y <- 1:(m*n)

  eta <- 2 * (1 + tanh(7 * (x - .8))) - 2
  eta0 <- eta
  for(i in 1:m){eta[group == i] <- eta[group == i] + avals[i]}
  for(i in 1:m){mu[group == i] <- 1 - 1 / (1 + exp(eta[group == i]))}
  for(i in 1:m){y[group == i] <- rbinom(n, size = 1, prob = mu[group == i])}
  dat <- data.frame(x = x, y = y, group = group)
  ansc <- cgamm(y ~ s.incr.conv(x) + (1|group),
  family = binomial(link = "logit"), reml = FALSE, data = dat)
  summary(ansc)
  anova(ansc)
}

\author{Xiyue Liao}

\keyword{cgamm routine}
