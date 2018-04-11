\name{cgam}

\alias{cgam}

\title{Constrained Generalized Additive Model Fitting}

\description{
  The partial linear generalized additive model is fitted using the method of maximum likelihood, where shape or order restrictions can be imposed on the non-parametrically modelled predictors with optional smoothing, and no restrictions are imposed on the optional parametrically modelled covariate. 
}

\usage{cgam(formula, nsim = 0, family = gaussian, cpar = 1.2, data = NULL, weights = NULL, 
sc_x = FALSE, sc_y = FALSE, pnt = TRUE, pen = 0)}

\arguments{
  \item{formula}{
    A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor". The response is a vector of length \eqn{n}. The specification of the model can be one of the three exponential families: gaussian, binomial and poisson. The systematic component \eqn{\eta} is \eqn{E(y)}, the log odds of \eqn{y = 1}, and the logarithm of \eqn{E(y)} respectively. A predictor can be a non-parametrically modelled variable with or without a shape or order restriction, or a parametrically modelled unconstrained covariate. In terms of a non-parametrically modelled predictor, the user is supposed to indicate the relationship between the systematic component \eqn{\eta} and a predictor \eqn{x} in the following way:

    Assume that \eqn{\eta} is the systematic component and \eqn{x} is a predictor: 
    \itemize{
    \item{incr(x): }{\eqn{\eta} is increasing in \eqn{x}. See \code{\link{incr}} for more details.}
    \item{s.incr(x): }{\eqn{\eta} is smoothly increasing in \eqn{x}. See \code{\link{s.incr}} for more details.}
    \item{decr(x): }{\eqn{\eta} is decreasing in \eqn{x}.  See \code{\link{decr}} for more details.}
    \item{s.decr(x): }{\eqn{\eta} is smoothly decreasing in \eqn{x}.  See \code{\link{s.decr}} for more details.}
    \item{conc(x): }{\eqn{\eta} is concave in \eqn{x}.  See \code{\link{conc}} for more details.}
    \item{s.conc(x): }{\eqn{\eta} is smoothly concave in \eqn{x}.  See \code{\link{s.conc}} for more details.}
    \item{conv(x): }{\eqn{\eta} is convex in \eqn{x}.  See \code{\link{conv}} for more details.}
    \item{s.conv(x): }{\eqn{\eta} is smoothly convex in \eqn{x}.  See \code{\link{s.conv}} for more details.}
    \item{incr.conc(x): }{\eqn{\eta} is increasing and concave in \eqn{x}.  See \code{\link{incr.conc}} for more details.}
    \item{s.incr.conc(x): }{\eqn{\eta} is smoothly increasing and concave in \eqn{x}. See \code{\link{s.incr.conc}} for more details.}
    \item{decr.conc(x): }{\eqn{\eta} is decreasing and concave in \eqn{x}. See \code{\link{decr.conc}} for more details.}
    \item{s.decr.conc(x): }{\eqn{\eta} is smoothly decreasing and concave in \eqn{x}. See \code{\link{s.decr.conc}} for more details.}
    \item{incr.conv(x): }{\eqn{\eta} is increasing and convex in \eqn{x}. See \code{\link{incr.conv}} for more details.}
    \item{s.incr.conv(x): }{\eqn{\eta} is smoothly increasing and convex in \eqn{x}. See \code{\link{s.incr.conv}} for more details.}
    \item{decr.conv(x): }{\eqn{\eta} is decreasing and convex in \eqn{x}. See \code{\link{decr.conv}} for more details.}
    \item{s.decr.conv(x): }{\eqn{\eta} is smoothly decreasing and convex in \eqn{x}. See \code{\link{s.decr.conv}} for more details.}
    \item{s(x): }{\eqn{\eta} is smooth in \eqn{x}. See \code{\link{s}} for more details.}
    \item{tree(x): }{\eqn{\eta} has a tree-ordering in \eqn{x}. See \code{\link{tree}} for more details.}
    \item{umbrella(x): }{\eqn{\eta} has an umbrella-ordering in \eqn{x}. See \code{\link{umbrella}} for more details.}
    }
}
  \item{nsim}{
  The number of simulations used to get the cic parameter. The default is nsim = 0. 
}
  \item{family}{
  A parameter indicating the error distribution and link function to be used in the model. It can be a character string naming a family function or the result of a call to a family function. This is borrowed from the glm routine in the stats package. There are three families used in cgam: gaussian, binomial, poisson, and Gamma. Note that if family = Ord is specified, a proportional odds regression model with shape constraints is fitted. This is under development. 
}

  \item{cpar}{A multiplier to estimate the model variance, which is defined as \eqn{\sigma^2 = SSR / (n - cpar * edf)}. SSR is the sum of squared residuals for the full model and edf is the effective degrees of freedom. The default is cpar = 1.2. The user-defined value must be between 1 and 2. See Meyer, M. C. and M. Woodroofe (2000) for more details.}

  \item{data}{
  An optional data frame, list or environment containing the variables in the model. The default is data = NULL.
  }
  \item{weights}{
  An optional non-negative vector of "replicate weights" which has the same length as the response vector. If weights are not given, all weights are taken to equal 1. The default is weights = NULL.
  }
  
  \item{sc_x}{Logical flag indicating if or not continuous predictors are normalized. The default is sc_x = FALSE.}
  \item{sc_y}{Logical flag indicating if or not the response variable is normalized. The default is sc_y = FALSE.}
  \item{pen}{User-defined penalty parameter. It must be non-negative. It will only be used in a warped-plane spline fit or a triangle spline fit. The default is pen = 0.}
  \item{pnt}{Logical flag indicating if or not penalized constrained regression splines are used. It will only be used in a warped-plane spline fit or a triangle spline fit. The default is pnt = TRUE.}
}

\details{
  We consider generalized partial linear models with independent observations from an exponential family of the form \eqn{p(y_i;\theta,\tau) = exp[\{y_i\theta_i - b(\theta_i)\}\tau - c(y_i, \tau)], i = 1,\ldots,n}, where the specifications of the functions \eqn{b} and \eqn{c} determine the sub-family of models. The mean vector \eqn{\mu = E(y)} has values \eqn{\mu_i = b'(\theta_i)}, and is related to a design matrix of predictor variables through a monotonically increasing link function \eqn{g(\mu_i) = \eta_i, i = 1,\ldots,n}, where \eqn{\eta} is the systematic component and describes the relationship with the predictors. The relationship between \eqn{\eta} and \eqn{\theta} is determined by the link function \eqn{b}.

For the additive model, the systematic component is specified for each observation by \eqn{\eta_i = f_1(x_{1i}) + \ldots + f_L(x_{Li}) + z_i'\beta}, where the functions \eqn{f_l} describe the relationships of the non-parametrically modelled predictors \eqn{x_l}, \eqn{\beta} is a parameter vector, and \eqn{z_i} contains the values of variables to be modelled parametrically. The non-parametric components are modelled with shape or order assumptions with optional smoothing, and the solution is obtained through an iteratively re-weighted cone projection, with no back-fitting of individual components.

Suppose that \eqn{\eta} is a \eqn{n} by \eqn{1} vector. The matrix form of the systematic component and the predictor is \eqn{\eta = \phi_1 + \ldots + \phi_L + Z\beta}, where \eqn{\phi_l} is the individual component for the \eqn{l}th non-parametrically modelled predictor, \eqn{l = 1, \ldots, L}, and \eqn{Z} is an \eqn{n} by \eqn{p} design matrix for the parametrically modelled covariate. 

To model the component \eqn{\phi_l}, smooth regression splines or non-smooth ordinal basis functions can be used. The constraints for the component \eqn{\phi_l} are in \eqn{C_l}. In the first case, \eqn{C_l} = \eqn{\{\phi_l \in R^n: \phi_l = v_l+B_l\beta_l}, where \eqn{\beta_l \ge 0} and \eqn{v_l\in V_l \}}, where \eqn{B_l} has regression splines as columns and \eqn{V_l} is the linear space contained in \eqn{C_l}, and in the second case, \eqn{C_l} = \eqn{\{\phi \in R^n: A_l\phi \ge 0} and \eqn{B_l\phi = 0\}}, for inequality constraint matrix \eqn{A_l} and equality constraint matrix \eqn{B_l}. 

The set \eqn{C_l} is a convex cone and the set \eqn{C = C_1 + \ldots + C_p + Z} is also a convex cone with a finite set of edges, where the edges are the generators of \eqn{C}, and \eqn{Z} is the column space of the design matrix \eqn{Z} for the parametrically modelled covariate. 

An iteratively re-weighted cone projection algorithm is used to fit the generalized regression model over the cone \eqn{C}. 

See references cited in this section and the official manual (\url{https://cran.r-project.org/package=coneproj}) for the R package coneproj for more details.
}

\value{
  %\item{vhat}{The fitted value of the systematic component \eqn{\eta} on the null space of the constraint set.}
  \item{etahat}{The fitted systematic component \eqn{\eta}.}
  \item{muhat}{The fitted mean value, obtained by transforming the systematic component \eqn{\eta} by the inverse of the link function.}
  \item{vcoefs}{The estimated coefficients for the basis spanning the null space of the constraint set.}
  \item{xcoefs}{The estimated coefficients for the edges corresponding to the smooth predictors with no shape constraint and shape-restricted predictors.}
  \item{zcoefs}{The estimated coefficients for the parametrically modelled covariate, i.e., the estimation for the vector \eqn{\beta}.}
  \item{ucoefs}{The estimated coefficients for the edges corresponding to the predictors with an umbrella-ordering constraint.}
  \item{tcoefs}{The estimated coefficients for the edges corresponding to the predictors with a tree-ordering constraint.}
  \item{coefs}{The estimated coefficients for the basis spanning the null space of the constraint set and edges corresponding to the shape-restricted and order-restricted predictors.}
  \item{cic}{The cone information criterion proposed in Meyer(2013a). It uses the "null expected degrees of freedom" as a measure of the complexity of the model. See Meyer(2013a) for further details of cic.}
  \item{d0}{The dimension of the linear space contained in the cone generated by all constraint conditions.}
  \item{edf0}{The estimated "null expected degrees of freedom". It is a measure of the complexity of the model. See Meyer (2013a) and Meyer (2013b) for further details.}
  \item{edf}{The constrained effective degrees of freedom.}
  \item{etacomps}{The fitted systematic component value for non-parametrically modelled predictors. It is a matrix of which each row is the fitted systematic component value for a non-parametrically modelled predictor. If there are more than one such predictors, the order of the rows is the same as the order that the user defines such predictors in the formula argument of cgam.}
  \item{y}{The response variable.}
  \item{xmat_add}{A matrix whose columns represent the shape-restricted predictors and smooth predictors with no shape constraint.}
  \item{zmat}{A matrix whose columns represent the basis for the parametrically modelled covariate. The user can choose to include a constant vector in it or not. It must have full column rank.}
  \item{ztb}{A list keeping track of the order of the parametrically modelled covariate.}
  \item{tr}{A matrix whose columns represent the predictors with a tree-ordering constraint.}
  \item{umb}{A matrix whose columns represent  the predictors with an umbrella-ordering constraint.}
  \item{tree.delta}{A matrix whose rows are the edges corresponding to the predictors with a tree-ordering constraint.}
  \item{umbrella.delta}{A matrix whose rows are the edges corresponding to the predictors with an umbrella-ordering constraint.}
  \item{bigmat}{A matrix whose rows are the basis spanning the null space of the constraint set and the edges corresponding to the shape-restricted and order-restricted predictors.}
  \item{shapes}{A vector including the shape and partial-ordering constraints in a cgam fit.}
  \item{shapesx}{A vector including the shape constraints in a cgam fit.}
 \item{prior.w}{User-defined weights.}
  \item{wt}{The weights in the final iteration of the iteratively re-weighted cone projections.}
 \item{wt.iter}{Logical flag indicating if or not iteratively re-weighted cone projections may be used. If the response is gaussian, then wt.iter = FALSE; if the response is binomial or poisson, then wt.iter = TRUE.}
  \item{family}{The family parameter defined in a cgam formula.}
  \item{SSE0}{The sum of squared residuals for the linear part.}
  \item{SSE1}{The sum of squared residuals for the full model.}
  \item{pvals.beta}{The approximate p-values for the estimation of the vector \eqn{\beta}. A t-distribution is used as the approximate distribution.}
  \item{se.beta}{The standard errors for the estimation of the vector \eqn{\beta}.}
  \item{null_df}{The degree of freedom for the null model of a cgam fit, i.e., the model only containing a constant vector.}
  \item{df}{The degree of freedom for the null space of a cgam fit.}
  \item{resid_df_obs}{The observed degree of freedom for the residuals of a cgam fit.}
  \item{null_deviance}{The deviance for the null model of a cgam fit, i.e., the model only containing a constant vector.}
  \item{deviance}{The residual deviance of a cgam fit.}
  \item{tms}{The terms objects extracted by the generic function \emph{terms} from a cgam fit. See the official help page (\url{http://stat.ethz.ch/R-manual/R-patched/library/stats/html/terms.html}) of the \emph{terms} function for more details.}
  \item{capm}{The number of edges corresponding to the shape-restricted predictors.}
  \item{capms}{The number of edges corresponding to the smooth predictors with no shape constraint.}
  \item{capk}{The number of non-constant columns of zmat.}
  \item{capt}{The number of edges corresponding to the tree-ordering predictors.}
  \item{capu}{The number of edges corresponding to the umbrella-ordering predictors.}
  \item{xid1}{A vector keeping track of the beginning position of the set of edges in bigmat for each shape-restricted predictor and smooth predictor with no shape constraint in xmat.}
  \item{xid2}{A vector keeping track of the end position of the set of edges in bigmat for each shape-restricted predictor and smooth predictor with no shape constraint in xmat.}
  \item{tid1}{A vector keeping track of the beginning position of the set of edges in bigmat for each tree-ordering factor in tr.}
  \item{tid2}{A vector keeping track of the end position of the set of edges in bigmat for each tree-ordering factor in tr.}
  \item{uid1}{A vector keeping track of the beginning position of the set of edges in bigmat for each umbrella-ordering factor in umb.}
  \item{uid2}{A vector keeping track of the end position of the set of edges in bigmat for each umbrella-ordering factor in umb.}
  \item{zid}{A vector keeping track of the positions of the parametrically modelled covariate.}
  \item{vals}{A vector storing the levels of each variable used as a factor.}
  \item{zid1}{A vector keeping track of the beginning position of the levels of each variable used as a factor.}
  \item{zid2}{A vector keeping track of the end position of the levels of each variable used as a factor.}
  \item{nsim}{The number of simulations used to get the cic parameter.}
  \item{xnms}{A vector storing the names of the shape-restricted predictors and the smooth predictors with no shape constraint in xmat.}
  \item{ynm}{The name of the response variable.}
  \item{znms}{A vector storing the names of the parametrically modelled covariate.}
  \item{is_param}{A logical scalar showing if or not a variable is a parametrically modelled covariate, which could be a linear term or a factor.}
  \item{is_fac}{A logical scalar showing if or not a variable is a factor.}
  \item{knots}{A list storing the knots used for each shape-restricted predictor and smooth predictor with no shape constraint. For a smooth,  constrained and a smooth, unconstrainted predictor, \emph{knots} is a vector of more than \eqn{1} elements, and for a shape-restricted predictor without smoothing, \emph{knots} = \eqn{0}.}
  \item{numknots}{A vector storing the number of knots for each shape-restricted predictor and smooth predictor with no shape constraint. For a smooth, constrained and a smooth, unconstrainted predictor, \emph{numknots} > \eqn{1}, and for a shape-restricted predictor without smoothing, \emph{numknots} = \eqn{0}.}
  \item{sps}{A character vector storing the \emph{space} parameter to create knots for each shape-restricted predictor.}
  \item{ms}{The centering terms used to make edges for shape-restricted predictors.}
\item{cpar}{The cpar argument in the cgam formula}
  \item{call}{The matched call.}
}

\references{
  Meyer, M. C. (2013a) Semi-parametric additive constrained regression. 
\emph{Journal of Nonparametric Statistics \bold{25(3)}}, 715

  Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications in Statistics \bold{42(5)}}, 1126--1139.  
  
  Meyer, M. C. and M. Woodroofe (2000) On the degrees of freedom in shape-restricted regression. \emph{Annals of Statistics \bold{28}}, 1083--1104.

  Meyer, M. C. (2008) Inference using shape-restricted regression splines. \emph{Annals of Applied Statistics \bold{2(3)}}, 1013--1033.

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
  fit <- cgam(Kyphosis ~ conc(Age) + incr.conc(Number) + decr.conc(Start), 
  family = binomial(), data = kyphosis) 
}

# Example 3.
  library(MASS)
  data(Rubber)
  
  # regress loss on hard and tens under the restrictions:
  # "decreasing" and "decreasing"
  fit.cgam <- cgam(loss ~ decr(hard) + decr(tens), data = Rubber)
  # "smooth and decreasing" and "smooth and decreasing"
  fit.cgam.s <- cgam(loss ~ s.decr(hard) + s.decr(tens), data = Rubber)
  summary(fit.cgam.s)
  anova(fit.cgam.s)
  
  # make a 3D plot based on fit.cgam and fit.cgam.s
  plotpersp(fit.cgam, th = 120, main = "3D Plot of a Cgam Fit")
  plotpersp(fit.cgam.s, tens, hard, data = Rubber, th = 120, main = "3D Plot of a Smooth Cgam Fit")
}

\keyword{cgam routine}

