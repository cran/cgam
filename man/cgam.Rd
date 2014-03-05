\name{cgam}

\alias{cgam}

\title{Constrained Generalized Additive Model Fitting}

\description{
  When the user provides a set of predictors with shape or order restriction, and a set of categorical covariates without any shape or order restriction, the cgam routine will return the maximum likelihood estimator for the constrained generalized additive model, using an iteratively re-weighted cone projection algorithm. The coneB routine in the R package coneproj is used in cgam to effect the hinge algorithm for one cone projection.
}

\usage{cgam(formula, nsim = 0, family = gaussian(), data = NULL, weights = NULL)}

\arguments{
  \item{formula}{
  A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor". The response is a vector of length \eqn{n}. The error distribution of the model can be one of the three exponential families: gaussian, binomial and poisson. A predictor can either be a non-parametrically modelled variable with a shape or order restriction, or a parametrically modelled unconstrained categorical covariate. The user can include both types of predictors in a cgam model. In terms of the shape or order restriction on a non-parametrically modelled predictor, the user is supposed to indicate the restriction by the following way:

    Assume that y is the response and x is a predictor: 
    \itemize{
    \item{incr(x): }{y is increasing in x.}
    \item{decr(x): }{y is decreasing in x.}
    \item{conc(x): }{y is concave in x.}
    \item{conv(x): }{y is convex in x.}
    \item{incr.conc(x): }{y is increasing and concave in x.}
    \item{decr.conc(x): }{y is decreasing and concave in x.}
    \item{incr.conv(x): }{y is increasing and convex in x.}
    \item{decr.conv(x): }{y is decreasing and convex in x.}
    \item{tree(x): }{y has a tree-ordering in x.}
    \item{umbrella(x): }{y has an umbrella-ordering in x.}
    }
}
  \item{nsim}{
  The number of simulations used to get the edf0 parameter. Note that when there is no shape-restricted or order-restricted predictor, nsim should be 0. The default is nsim = 0.
}
  \item{family}{
  A parameter indicating the error distribution and link function to be used in the model. It can be a character string naming a family function or the result of a call to a family function. This is borrowed from the glm routine in the stats package. There are three families used in cgam: gaussian, binomial and poisson.
}
  \item{data}{
  An optional data frame, list or environment containing the variables in the model. The default is data = NULL.
  }
  \item{weights}{
  An optional non-negative vector of "prior weights" which has the same length as the response vector. If weights is not given, all weights are taken to equal 1. The default is weights = NULL.
  }

}

\details{
  We consider generalized partial linear models with independent observations from an exponential family of the form \eqn{p(y_i;\theta,\tau) = exp[\{y_i\theta_i - b(\theta_i)\}\tau - c(y_i, \tau)], i = 1,\ldots,n}, where the specifications of the functions \eqn{b} and \eqn{c} determine the sub-family of models. The mean vector \eqn{\mu = E(y)} has values \eqn{\mu_i = b'(\theta_i)}, and is related to a design matrix of predictor variables through a monotonically increasing link function \eqn{g(\mu_i) = \eta_i, i = 1,\ldots,n}, where \eqn{\eta} is the "systematic component" and 
describes the relationship with the predictors. The relationship between \eqn{\eta} and \eqn{\theta} is determined by the link function \eqn{b}.

For the additive model, the systematic component is specified for each observation by \eqn{\eta_i = f_1(x_{1i}) + \ldots + f_L(x_{Li}) + \beta ' z_i}, where the functions \eqn{f_l} describe the relationships of the non-parametrically modelled predictors \eqn{x_l}, \eqn{\beta} is a parameter vector, and \eqn{z_i} contains the values of variables to be modelled parametrically. The non-parametric components are modelled with shape or order assumptions, and the solution is obtained through an iteratively re-weighted cone projection, with no back-fitting of individual components.

Suppose that \eqn{\eta} is a \eqn{n} by \eqn{1} vector. The matrix form of the systematic component and the predictor is \eqn{\eta = \phi_1 + \ldots + \phi_L + Z\beta}, where \eqn{\phi_l} is the individual component for the \eqn{l}th non-parametrically modelled predictor, \eqn{l = 1, \ldots, L}, and \eqn{Z} is an \eqn{n} by \eqn{p} design matrix for the parametically modelled predictors. 

The constraints for the component \eqn{\phi_l} are in \eqn{C_l}, where \eqn{C_l} = \eqn{\{\phi: A_l\phi \ge 0} and \eqn{B_l\phi = 0\}}, for matrices \eqn{A_l} and \eqn{B_l}. 

The set \eqn{C_l} is a convex cone and the set \eqn{C = C_1 + \ldots + C_p + Z} is also a convex cone with a finite set of edges, where the edges are the generator of \eqn{C} and \eqn{Z} is the column space of the design matrix \eqn{Z} for the parametrically modelled predictors. 

An iteratively re-weighted cone projection algorithm is used to fit the generalized regression model over the cone \eqn{C}. 

See references cited in this section and the official manual for the R package coneproj for more details.
}

\value{
  \item{vhat}{The fitted value of the response on the null space of the constraint set.}
  \item{etahat}{The fitted systematic component \eqn{\eta}.}
  \item{muhat}{The fitted mean value, obtained by transforming the systematic component \eqn{\eta} by the inverse of the link function.}
  \item{vcoefs}{The estimated coefficients for the basis spanning the null space of the constraint set.}
  \item{xcoefs}{The estimated coefficients for the edges corresponding to the shape-restricted predictors.}
  \item{zcoefs}{The estimated coefficients for zmat, i.e., the estimation for the vector \eqn{\beta}.}
  \item{ucoefs}{The estimated coefficients for the edges corresponding to the predictors with an umbrella-ordering constraint.}
  \item{tcoefs}{The estimated coefficients for the edges corresponding to the predictors with a tree-ordering constraint.}
  \item{coefs}{The estimated coefficients for the basis spanning the null space of the constraint set and edges corresponding to the constrained predictors.}
  \item{cic}{The cone information criterion proposed in Meyer(2013a). It uses the "null expected degrees of freedom" as a measure of the complexity of the model. See Meyer(2013a) for further details of cic.}
  \item{d0}{The dimension of the null space contained in the cone generated by all constraint conditions.}
  \item{edf0}{The estimated "null expected degrees of freedom". It is a measure of the complexity of the model. See Meyer (2013a) and Meyer (2013b) for further details.}
  \item{etacomps}{The fitted mean value for the shape-restricted or order-restricted predictors. It is a matrix of which each row is the fitted mean value for a shape-restricted or order-restricted predictor. If there are more than one such predictor, the order of the rows is the same as the order that the user defines shape-restricted or order-restricted predictors in the formula argument of cgam.}
  \item{xmat}{A matrix whose columns respresent the shape-restricted predictors.}
  \item{zmat}{A matrix whose columns represent the unconstrained categorical covariates. The user can choose to include a constant vector in it or not. It must be of full column rank.}
  \item{tr}{A matrix whose columns respresent the predictors with a tree-ordering constraint.}
  \item{umb}{A matrix whose columns respresent  the predictors with an umbrella-ordering constraint.}
  \item{tree.delta}{A matrix whose rows are the edges corresponding to the predictors with a tree-ordering constraint.}
  \item{umbrella.delta}{A matrix whose rows are the edges corresponding to the predictors with an umbrella-ordering constraint.}
  \item{bigmat}{A matrix whose rows are the edges corresponding to the constrained predictors.}
  \item{shapes}{A vector including the shape constraints in a cgam fit.}
  \item{wt}{The weights in the final iteration of the iteratively re-weighted cone projections.}
 \item{wt.iter}{A logical scalar indicating if or not iteratively re-weighted cone projections are used to get the fit. If the error distribution is gaussian, then wt.iter = FALSE; if the error distribution is binomial or poisson, then wt.iter = TRUE.}
  \item{family}{The family parameter defined by the user in cgam. }
  \item{SSE0}{The sum of squared residuals for the linear part.}
  \item{SSE1}{The sum of squared residuals for the full model.}
  \item{pvals.beta}{The approximate p-values for the estimation of the vector \eqn{\beta}. A t-distribution is used as the approximate distribution.}
  \item{se.beta}{The standard errors for the estimation of the vector \eqn{\beta}.}
  \item{null_df}{The degree of freedom for the null model of a cgam fit, i.e., the model only containing a constant vector.}
  \item{df}{The degree of freedom for the null space of a cgam fit.}
  \item{resid_df_obs}{The observed degree of freedom for the residuals of a cgam fit.}
  \item{null_deviance}{The deviance for the null model of a cgam fit, i.e., the model only containing a constant vector.}
  \item{deviance}{The residual deviance of a cgam fit.}
  \item{tms}{The terms objects extracted by the generic function \emph{terms} from a cgam fit. See the help page of the \emph{terms} function for more details.}
  \item{capm}{The number of edges corresponding to the shape-restricted predictors.}
  \item{capk}{The number of non-constant columns of zmat.}
  \item{capt}{The number of edges corresponding to the tree-ordering predictors.}
  \item{capu}{The number of edges corresponding to the umbrella-ordering predictors.}
  \item{xid1}{A vector storing the indices keeping track of the beginning position of the set of edges for each shape-restricted predictor in xmat.}
  \item{xid2}{A vector storing the indices keeping track of the end position of the set of edges for each shape-restricted predictor in xmat.}
  \item{tid1}{A vector storing the indices keeping track of the beginning position of the set of edges for each tree-ordering factor in tr.}
  \item{tid2}{A vector storing the indices keeping track of the end position of the set of edges for each tree-ordering factor in tr.}
  \item{uid1}{A vector storing the indices keeping track of the beginning position of the set of edges for each umbrella-ordering factor in umb.}
  \item{uid2}{A vector storing the indices keeping track of the end position of the set of edges for each umbrella-ordering factor in umb.}
  \item{call}{The matched call.}
}

\references{
  Meyer, M. C. (2013a) Semi-parametric additive constrained regression. 
\emph{Journal of Nonparametric Statistics \bold{25(3)}}, 715

  Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications in Statistics \bold{42(5)}}, 1126--1139.  
  
  Meyer, M. C. and M. Woodroofe (2000) On the degrees of freedom in shape-restricted regression. \emph{Annals of Statistics \bold{28}}, 1083--1104.

  Mammen, E. and K. Yu (2007) Additive isotonic regression. \emph{IMS Lecture  Notes-Monograph Series Asymptotics: Particles, Process, and Inverse Problems \bold{55}}, 179--195.

  Huang, J. (2002) A note on estimating a partly linear model under monotonicity constraints. \emph{Journal of Statistical Planning and Inference \bold{107}}, 343--351.

  Cheng, G.(2009) Semiparametric additive isotonic regression. \emph{Journal of Statistical Planning and Inference \bold{139}}, 1980--1991.

  Bacchetti, P. (1989) Additive isotonic models. \emph{Journal of the American Statistical Association \bold{84(405)}}, 289--294.

}

\author{Mary C. Meyer and Xiyue Liao }

\examples{
# Example 1.
  data(cubic)
  # extract x
  x <- cubic$x

  # extract y
  y <- cubic$y

  # regress y on x with no restriction with lm()
  fit.lm <- lm(y ~ x + I(x^2) + I(x^3))

  # regress y on x under the restriction: "increasing and convex"
  fit.cgam <- cgam(y ~ incr.conv(x))

  # make a plot to compare the two fits
  par(mar = c(4, 4, 1, 1))
  plot(x, y, cex = .7, xlab = "x", ylab = "y")
  lines(x, fit.cgam$muhat, col = 2, lty = 2)
  lines(x, fitted(fit.lm), col = 1, lty = 1)
  legend("topleft", bty = "n", c("constrained cgam fit", "unconstrained lm fit"), 
  lty = c(2, 1), col = c(2, 1))

# Example 2.
\dontrun{
  library(gam)
  data(kyphosis)
  
  # regress Kyphosis on Age, Number, and Start under the restrictions:
  # "concave", "increasing and concave", and "decreasing and concave" 
  fit <- cgam(Kyphosis ~ conc(Age) + incr.conc(Number) + decr.conc(Start), family = binomial(), 
  nsim = 1e+3, data = kyphosis) 
}

# Exampe 3.
  library(MASS)
  data(Rubber)
  
  # regress loss on hard and tens under the restrictions:
  # "decreasing" and "decreasing"
  fit.cgam <- cgam(loss ~ decr(hard) + decr(tens), data = Rubber)
  n <- 10
  hard <-  Rubber$hard
  tens <- Rubber$tens

  # make a 3D plot based on fit.cgam 
  x1grid <- seq(min(hard), max(hard), length = n)
  x2grid <- seq(min(tens), max(tens), length = n)	
  TwoDdata <- data.frame(hard = rep(x1grid, each = length(x2grid)), 
  tens = rep(x2grid, length(x1grid)))
  thint <- predict(fit.cgam, TwoDdata)$fit
  A <- matrix(thint, length(x1grid), length(x2grid), byrow = TRUE) 
  persp(A, xlab = "hard", ylab = "tens", zlab = "loss", theta = 120)
  title("3D Plot of a Cgam Fit")
}

\keyword{cgam routine}

