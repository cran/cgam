\name{wps}

\alias{wps}

\title{Warped-Plane Spline-Based Isotonic Regression}

\description{
Isotonic regression surface in two dimensions using warped-plane splines is fitted by the wps routine without using additivity assumptions. Meanwhile, linear covariate effects can also be incorporated in the model. The surface and covariate effects are estimated simultaneously with a single cone projection. A penalized spline fit is also provided which allows for more knots and greater flexibility.

}

\usage{
wps(formula, family = gaussian(), data = NULL, weights = NULL, pnt = FALSE, 
pen = 0, cpar = 1.5)
}

\arguments{
  \item{formula}{ A formula object which gives a symbolic description of the model to be fitted. It has the form "response ~ predictor". The response is a vector of length \eqn{n}. The user is supposed to indicate the relationship between the mean value \eqn{E(y)} and the two predictors \eqn{x_1} and \eqn{x_2} in the following way:

    Assume that \eqn{E(y)} is the mean value and we model the relationship between \eqn{E(y)} and the two predictors \eqn{x_1} and \eqn{x_2} as: 
    \itemize{
    \item{dd(x1, x2): }{\eqn{E(y)} is decreasing in both \eqn{x_1} and \eqn{x_2}. See \code{\link{dd}} for more details.}
    \item{di(x1, x2):  }{\eqn{E(y)} is decreasing in \eqn{x_1} and increasing in \eqn{x_2}. See \code{\link{di}} for more details.}
    \item{ii(x1, x2): }{\eqn{E(y)} is increasing in both \eqn{x_1} and \eqn{x_2}.  See \code{\link{ii}} for more details.}
    }
}

  \item{family}{A parameter indicating the error distribution and link function to be used in the model. It can be a character string naming a family function or the result of a call to a family function. This is borrowed from the glm routine in the stats package. For now, only the option family = gaussian() is allowed in the wps routine.}

  \item{data}{An optional data frame, list or environment containing the variables in the model. The default is data = NULL.}

  \item{weights}{An optional non-negative vector of "replicate weights" which has the same length as the response vector. If weights are not given, all weights are taken to equal 1. The default is weights = NULL.}

  \item{pnt}{Logical flag indicating if or not penalized warped-plane splines are used. When pnt is TRUE, penalized warped-plane splines are used; otherwise, unpenalized warped-plane splines are used. The default is pnt = FALSE. Note that if pnt is TRUE, the user must define a positive value for the argument pen.}

  \item{pen}{The penalty term for penalized warped-plane splines. It can only be a non-negative real number. When the user chooses penalized warped-plane splines, a positive penalty term must be defined. If this term is zero, then unpenalized warped-plane splines are used. The default is pen = \eqn{0}.}

  \item{cpar}{A multiplier to estimate the model variance, which is defined as \eqn{\sigma^2 = SSE / (n - cpar * edf)}. SSE is the sum of squared residuals for the full model and edf is the effective degrees of freedom. The default is cpar = 1.5. See Meyer, M. C. and M. Woodroofe (2000) for more details.}
}

\details{
  We consider the regression model  \eqn{y_i = f(t_{1i}, t_{2i}) + z_i'\beta +\varepsilon_{i}, i = 1,\ldots,n}, where \eqn{\beta} is a \eqn{p} by \eqn{1} parameter vector, and the \eqn{\varepsilon_i's} are mean-zero random errors. We know \emph{a priori} that \eqn{f} is continuous and isotonic in both dimensions; that is, for fixed \eqn{t_{1}} and \eqn{z} values, if \eqn{t_{2a} \leq t_{2b}}, then \eqn{f(t_{1}, t_{2a}) \leq f(t_{1}, t_{2a})}, and similarly for fixed \eqn{t_{2}} and \eqn{z} values, \eqn{f} is non-decreasing in \eqn{t_{1}}. For splines of degree two or higher, obtaining a finite set of linear inequality constraints that are necessary and sufficient for isotonicity in both dimensions does not seem to be feasible. However, if we use linear spline basis functions, then the necessary and sufficient constraints are straight-forward and the fitted surface can be described as a continuous piece-wise warped plane, called a "warped-plane spline" (WPS) surface. 

The surface and covariate effects are estimated simultaneously with a single cone projection (no back-fitting). See references cited in this section and the official manual(\url{https://cran.r-project.org/package=coneproj})for the R package coneproj for more details.

A penalized spline version is also provided. Over each knot rectangle, the regression surface is a warped plane, and the slopes can change abruptly from one rectangle to the next. To obtain smoother fits, and to side-step the problem of knot choices, we can use large number of knots for both predictors and penalize these changes in slopes. The size of the penalty parameter will control the effective degrees of freedom of the fit. In practice, the penalty term can be chosen through generalized cross-validation, similar to the method in Meyer (2012).


}

\value{
\item{k1}{Knots used for \eqn{x_1}.}
\item{k2}{Knots used for \eqn{x_2}.}
\item{muhat}{The estimated constrained mean value.}
\item{muhatu}{The estimated unconstrained mean value.}
\item{SSE1}{The sum of squared residuals for the full model.}
\item{SSE0}{The sum of squared residuals for the linear part.}
\item{edf}{The constrained effective degrees of freedom.}
\item{edfu}{The unconstrained effective degrees of freedom.}
\item{delta}{A matrix whose columns are the edges corresponding to the isotonically modelled predictors \eqn{x_1} and \eqn{x_2}.}
\item{zmat}{A matrix whose columns represent the parametrically modelled covariate. The user can choose to include a constant vector in it or not. It must have full column rank.}
\item{xmat}{A matrix whose columns are \eqn{x_1} and \eqn{x_2}.}
\item{coefs}{The estimated constrained coefficients for the basis spanning the null space of the constraint set and edges corresponding to isotonically modelled predictors \eqn{x_1} and \eqn{x_2}.}
\item{coefsu}{The estimated unconstrained coefficients for the basis spanning the null space of the constraint set and edges corresponding to \eqn{x_1} and \eqn{x_2}.}
\item{zcoefs}{The estimated coefficients for the parametrically modelled covariate.}
%\item{sighat}{The estimated model variance.}
\item{pvals.beta}{The approximate p-values for the estimation of the vector \eqn{\beta}. A t-distribution is used as the approximate distribution.}
\item{se.beta}{The standard errors for the estimation of the vector \eqn{\beta}.}
%\item{covmat}{The estimated covariance matrix of \eqn{\hat{\beta}}.}
\item{gcv}{The generalized cross validation (GCV) value for the constrained fit.}
\item{gcvu}{The generalized cross validation (GCV) value for the unconstrained fit.}
\item{xnms}{A vector storing the names of \eqn{x_1} and \eqn{x_2}.}
\item{znms}{A vector storing the names of the parametrically modelled covariate.}
\item{zid}{A vector keeping track of the positions of the parametrically modelled covariate.}
\item{vals}{A vector storing the levels of each variable used as a factor.}
\item{zid1}{A vector keeping track of the beginning position of the levels of each variable used as a factor.}
\item{zid2}{A vector keeping track of the end position of the levels of each variable used as a factor.}
\item{ynm}{The name of the response variable.}
\item{decrs}{A vector of two logical values indicating the monotonicity of the isotonically-constrained surface with respect to \eqn{x_1} and \eqn{x_2}.}
\item{tms}{The terms objects extracted by the generic function \emph{terms} from a wps fit. See the official help page (\url{http://stat.ethz.ch/R-manual/R-patched/library/stats/html/terms.html}) of the \emph{terms} function for more details.}
\item{is_param}{A logical scalar showing if or not a variable is a parametrically modelled covariate, which could be a linear term or a factor.}
\item{is_fac}{A logical scalar showing if or not a variable is a factor.}
\item{call}{The matched call.}
}

\author{
  Mary C. Meyer and Xiyue Liao
}

\references{
Meyer, M. C. and M. Woodroofe (2000) On the degrees of freedom in shape-restricted regression. \emph{Annals of Statistics \bold{28}}, 1083--1104.

Meyer, M. C. (2012) Constrained penalized splines. \emph{Canadian Journal of Statistics \bold{40(1)}}, 190--206.

Meyer, M. C. (2016) Estimation and inference for isotonic regression in two dimensions, using warped-plane splines.

}


\examples{
  library(MASS)
  data(Rubber)

  # regress loss on hard and tens under the shape-restriction: "doubly-decreasing" 
  # with a penalty term equal to 1
  # use 13 knots for each predictor
  ans <- wps(loss ~ dd(hard, tens, numknots = c(13, 13)), data = Rubber, pen = 1)

  # make a 3D plot of the constrained surface
  plotpersp(ans, hard, tens, data = Rubber)
}

\keyword{wps routine}

