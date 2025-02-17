\name{ShapeSelect}
\alias{ShapeSelect}

\title{Variable and Shape Selection via Genetic Algorithm}

\description{
The partial linear generalized additive model is considered, where the goal is to choose a subset of predictor variables and describe the component relationships with the response, in the case where there is very little \emph{a priori} information. For each predictor, the user need only specify a set of possible shape or order restrictions. A model selection method chooses the shapes and orderings of the relationships as well as the variables.  For each possible combination of shapes and orders for the predictors, the maximum likelihood estimator for the constrained generalized additive model is found using iteratively re-weighted cone projections.  The cone information criterion is used to select the best combination of variables and shapes.  
}

\usage{ShapeSelect(formula, family = gaussian, cpar = 2, data = NULL, weights = NULL, 
npop = 200, per.mutate = 0.05, genetic = FALSE)}

\arguments{
  \item{formula}{

 A formula object which includes a set of predictors to be selected. It has the form "response ~ predictor". The response is a vector of length \eqn{n}. The specification of the model can be one of the three exponential families: gaussian, binomial and poisson. The systematic component \eqn{\eta} is \eqn{E(y)}, the log odds of \eqn{y = 1}, and the logarithm of \eqn{E(y)} respectively. The user is supposed to define at least one predictor in the formula, which could be a non-parametrically modelled variable or a parametrically modelled covariate (categorical or linear). Assume that \eqn{\eta} is the systematic component and \eqn{x} is a predictor, two symbolic routines \code{\link{shapes}} and \code{\link{in.or.out}} are used to include \eqn{x} in the formula. 
    \describe{
    \item{shapes(x): }{\eqn{x} is included as a non-parametrically modelled predictor in the formula. See \code{\link{shapes}} for more details.}
    \item{in.or.out(x): }{\eqn{x} is included as a categorical or linear predictor in the formula. See \code{\link{in.or.out}} for more details.}
    }
  }
  \item{family}{A parameter indicating the error distribution and link function to be used in the model. It can be a character string naming a family function or the result of a call to a family function. This is borrowed from the glm routine in the stats package. There are three families used in ShapeSelect: gaussian, binomial and poisson.}
  \item{cpar}{A multiplier to estimate the model variance, which is defined as \eqn{\sigma^2 = SSR / (n - d_0 - cpar * edf)}. SSR is the sum of squared residuals for the full model, \eqn{d_0} is the dimension of the linear space in the cone, and edf is the effective degrees of freedom. The default is cpar = 2. See Meyer, M. C. and M. Woodroofe (2000) for more details.}
  \item{data}{An optional data frame, list or environment containing the variables in the model. The default is data = NULL.}
  \item{weights}{An optional non-negative vector of "replicate weights" which has the same length as the response vector. If weights are not given, all weights are taken to equal 1. The default is weights = NULL.}
  \item{npop}{The population size used for the genetic algorithm. The default is npop = 200.}
  \item{per.mutate}{The percentage of mutation used for the genetic algorithm. The default is per.mutate = 0.05}
  \item{genetic}{A logical scalar showing if or not the genetic algorithm is defined by the user to select the best model. The default is genetic = FALSE.}
}

\details{

Note that when the argument genetic is set to be FALSE, the routine will check to see if using the genetic algorithm is better than going through all models to find the best fit. The primary concern is running time. An interactive dialogue window may pop out to ask the user if they prefer to the genetic algorithm when it may take too long to do a brutal search, and if there are too many possible models to go through, like more than one million, the routine will implement the genetic algorithm anyway.

See references cited in this section for more details.
}

\value{
\item{top}{The best model of the final population, which shows the variables chosen along with their best shapes.}
\item{pop}{The final population ordered by its fitness values. It is a data frame, and each row of this data frame shows the shapes chosen for predictors in an individual model. Besides, the fitness value for each individual model is included in the last column of the data frame. For example, we have two continuous predictors \eqn{x_1}, \eqn{x_2}, and a categorical predictor \eqn{z}, then a row of this data frame may look like: "flat", "s.incr", "in", \eqn{-12.3806}, which means that \eqn{x_1} is not chosen, \eqn{x_2} is chosen with the shape constraint to be smoothly increasing, \eqn{z} is included in the model, and the fitness value for the model is \eqn{-12.3806}.}
\item{fitness}{The sorted fitness values for the final population.}
\item{tm}{Total cpu running time.}
\item{xnms}{A vector storing the name of the nonparametrically-modelled predictor in a ShapeSelect formula.}
\item{znms}{A vector storing the name of the parametrically-modelled predictor in a ShapeSelect formula, which is a categorical predictor or a linear term.}
\item{trnms}{A vector storing the name of the treatment predictor in a ShapeSelect formula, which has three possible levels: no effect, tree ordering, unordered.}
\item{zfacs}{A logical vector keeping track of if the parametrically-modelled predictor in a ShapeSelect formula is a categorical predictor or a linear term.}
\item{mxf}{A vector keeping track of the largest fitness value in each generation.}
\item{mnf}{A vector keeping track of the mean fitness value in each generation.}
\item{GA}{A logical scalar showing if or not the genetic algorithm is actually implemented to select the best model.}
\item{best.fit}{The best model fitted by the cgam routine, given the best variables with their shape constraints chosen by the ShapeSelect routine.}
\item{call}{The matched call.}
}

\references{
 Meyer, M. C. (2013a) Semi-parametric additive constrained regression. \emph{Journal of Nonparametric Statistics \bold{25(3)}}, 715

  Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications in Statistics \bold{42(5)}}, 1126--1139.  
  
  Meyer, M. C. and M. Woodroofe (2000) On the degrees of freedom in shape-restricted regression. \emph{Annals of Statistics \bold{28}}, 1083--1104.

  Meyer, M. C. (2008) Inference using shape-restricted regression splines. \emph{Annals of Applied Statistics \bold{2(3)}}, 1013--1033.

  Meyer, M. C. and Liao, X. (2016) Variable and shape selection for the generalized additive model. \emph{under preparation}
}

\author{Mary C. Meyer and Xiyue Liao}

\seealso{\code{\link{cgam}}}

\examples{
\dontrun{
# Example 1.
  library(MASS)
  data(Rubber)

  # ShapeSelect can be used to go through all models to find the best model
  fit <- ShapeSelect(loss ~ shapes(hard, set = "s.9") + shapes(tens, set = "s.9"),
  data = Rubber, genetic = FALSE)
 
  # the user can also choose to find the best model by the genetic algorithm
  # given any total number of possible models
  fit <- ShapeSelect(loss ~ shapes(hard, set = "s.9") + shapes(tens, set = "s.9"),
  data = Rubber, genetic = TRUE)

  # check the best model
  fit$top

  # check the running time
  fit$tm

# Example 2.
  # simulate a data set such that the mean is smoothly increasing-convex in x1 and x2
  n <- 100
  x1 <- runif(n)
  x2 <- runif(n)
  y0 <-  x1^2 + x2 + x2^3

  z <- rep(0:1, 50)
  for (i in 1:n) {
    if (z[i] == 1) 
      y0[i] = y0[i] * 1.5
  }

  # add some random errors and call the routine
  y <- y0 + rnorm(n)

  # include factor(z) in the formula and determine if factor(z) should be in the model 
  fit <- ShapeSelect(y ~ shapes(x1, set = "s.9") + shapes(x2, set = "s.9") + in.or.out(factor(z)))
  
  # use the genetic algorithm
  fit <- ShapeSelect(y ~ shapes(x1, set = "s.9") + shapes(x2, set = "s.9") + in.or.out(factor(z)),
   npop = 300, per.mutate=0.02)

  # include z as a linear term in the formula and determine if z should be in the model 
  fit <- ShapeSelect(y ~ shapes(x1, set = "s.9") + shapes(x2, set = "s.9") + in.or.out(z))

  # include z as a linear term in the model 
  fit <- ShapeSelect(y ~ shapes(x1, set = "s.9") + shapes(x2, set = "s.9") + z)

  # include factor(z) in the model 
  fit <- ShapeSelect(y ~ shapes(x1, set = "s.9") + shapes(x2, set = "s.9") + factor(z))

  # check the best model
  bf <- fit$best.fit 
 
  # make a 3D plot of the best fit
  plotpersp(bf, categ = "z")
}
}

\keyword{variable and shape selection}
\keyword{genetic algorithm}
