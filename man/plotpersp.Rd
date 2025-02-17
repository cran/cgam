\name{plotpersp}

\alias{plotpersp}

\title{Create a 3D Plot for a CGAM Object}

\description{Given an object of the cgam class, which has at least two non-parametrically modelled predictors, this routine will make a 3D plot of the fit with a set of two non-parametrically modelled predictors in the formula being the x and y labs. If there are more than two non-parametrically modelled predictors, any other such predictor will be evaluated at the largest value which is smaller than or equal to its median value.

If there is any categorical covariate and if the user specifies the argument categ to be a character representing a categorical covariate in the formula, then a 3D plot with multiple parallel surfaces, which represent the levels of a categorical covariate in an ascending order, will be created; otherwise, a 3D plot with only one surface will be created. Each level of a categorical covariate will be evaluated at its mode.

This routine is extended to make a 3D plot for an object fitted by warped-plane splines or triangle splines. Note that two non-parametrically modelled predictors specified in this routine must both be modelled as addtive components, or a pair of predictors forming an isotonic or convex surface without additivity assumption.

This routine is an extension of the generic R graphics routine persp. See the documentation below for more details. 

}

\usage{
  plotpersp(object, x1 = NULL, x2 = NULL, ...)
}

\arguments{
  \item{object}{An object of the cgam class with at least two non-parametrically modelled predictors.}
  \item{x1}{A non-parametrically modelled predictor in a cgam fit. If the user omits x1 and x2, then the first two non-parametric predictors in a cgam formula will be used.}
  \item{x2}{A non-parametrically modelled predictor in a cgam fit. If the user omits x1 and x2, then the first two non-parametric predictors in a cgam formula will be used.}
\item{...}{
Arguments to be passed to the S3 method for the cgam class:
\describe{
   \item{x1nm: }{Character name of x1.}
  \item{x2nm: }{Character name of x2.}
  \item{data: }{The data frame based on which the user get a cgam fit.}
  \item{surface: }{The type of the surface of a 3D plot. For a cgam fit, if surface == "mu", then the surface of the estimated mean value of the fit will be plotted; if surface == "eta", then the surface of the estimated systematic component value of the fit will be plotted. The default is surface = "mu"; for a warped-plane spline fit, if surface == "C", then the surface of the constrained estimated mean value of the fit will be plotted, while if surface == "U", then the surface of the unconstrained estimated mean value of the fit will be plotted. The default is surface = "C".}
 \item{categ: }{Optional categorical covariate(s) in a cgam fit. If there is any categorical covariate and if the user specifies the argument categ to be a character representing a categorical covariate in the formula, then a 3D plot with multiple parallel surfaces, which represent the levels of a categorical covariate in an ascending order, will be created; otherwise, a 3D plot with only one surface will be created. Each level of a categorical covariate will be evaluated at its mode. The default is categ = NULL.}
  \item{col: }{The color(s) of a 3D plot created by plotpersp. If col == NULL, "white" will be used when there is only one surface in the plot, and a sequence of colors will be used in a fixed order when there are multiple parallel surfaces in the plot. For example, when there are two surfaces, the lower surface will be in the color "peachpuff", and the higher surface will be in the color "lightblue". The default is col =  NULL.}
  \item{random: }{A logical scalar. If random == TRUE, color(s) for a 3D plot will be randomly chosen from ten colors, namely, "peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke"; otherwise, "white" will be used when there is only one surface in the plot, and a sequence of colors will be used in a fixed order when there are multiple parallel surfaces in the plot.}
  \item{ngrid: }{This is a positive integer specifying how dense the \eqn{x} grid and the \eqn{y} grid will be. The default is ngrid = \eqn{12}. Note that this argument is only used for a cgam fit.}
  \item{xlim: }{The xlim argument inherited from the persp routine.}
  \item{ylim: }{The ylim argument inherited from the persp routine.}
  \item{zlim: }{The zlim argument inherited from the persp routine.}
  \item{xlab: }{The xlab argument inherited from the persp routine.}
  \item{ylab: }{The ylab argument inherited from the persp routine.}
  \item{zlab: }{The zlab argument inherited from the persp routine.}
  \item{main: }{The main argument inherited from the persp routine.}
  \item{th: }{The theta argument inherited from the persp routine.}
  \item{ltheta: }{The ltheta argument inherited from the persp routine.}
  \item{main: }{The main argument inherited from the persp routine.}
  \item{ticktype: }{The ticktype argument inherited from the persp routine.}
  %\item{\dots}{Further arguments passed to the routine.}
}
}
}

\value{The routine plotpersp returns a 3D plot of an object of the cgam class. The \eqn{x} lab and \eqn{y} lab represent a set of non-parametrically modelled predictors used in a cgam formula, and the \eqn{z} lab represents the estimated mean value or the estimated systematic component value.}

\author{
 Mary C. Meyer and Xiyue Liao
}

\examples{
# Example 1.
  data(FEV)

  # extract the variables
  y <- FEV$FEV
  age <- FEV$age
  height <- FEV$height
  sex <- FEV$sex
  smoke <- FEV$smoke

  fit <- cgam(y ~ incr(age) + incr(height) + factor(sex) + factor(smoke), nsim = 0)
  fit.s <- cgam(y ~ s.incr(age) + s.incr(height) + factor(sex) + factor(smoke), nsim = 0)

  plotpersp(fit, age, height, ngrid = 10, main = "Cgam Increasing Fit", 
  sub = "Categorical Variable: Sex", categ = "factor(sex)")
  plotpersp(fit.s, age, height, ngrid = 10, main = "Cgam Smooth Increasing Fit", 
  sub = "Categorical Variable: Smoke", categ = "factor(smoke)")

# Example 2.
  data(plasma)

  # extract the variables
  y <- plasma$logplasma
  bmi <- plasma$bmi
  logdietfat <- plasma$logdietfat
  cholest <- plasma$cholest
  fiber <- plasma$fiber
  betacaro <- plasma$betacaro
  retinol <- plasma$retinol
  smoke <- plasma$smoke
  vituse <- plasma$vituse

  fit <- cgam(y ~  s.decr(bmi) + s.decr(logdietfat) + s.decr(cholest) + s.incr(fiber) 
+ s.incr(betacaro) + s.incr(retinol) + factor(smoke) + factor(vituse)) 

  plotpersp(fit, bmi, logdietfat, ngrid = 15, th = 120, ylab = "log(dietfat)", 
zlab = "est mean of log(plasma)", main = "Cgam Fit with the Plasma Data Set", 
sub = "Categorical Variable: Vitamin Use", categ = "factor(vituse)")

# Example 3.
  data(plasma)
  addl <- 1:314*0 + 1 
  addl[runif(314) < .3] <- 2
  addl[runif(314) > .8] <- 4
  addl[runif(314) > .8] <- 3

  ans <- cgam(logplasma ~ s.incr(betacaro, 5) + s.decr(bmi) + s.decr(logdietfat) 
+ as.factor(addl), data = plasma)
  plotpersp(ans, betacaro, logdietfat, th = 240, random = TRUE, 
categ = "as.factor(addl)", data = plasma)

# Example 4.
\dontrun{
  n <- 100
  set.seed(123)
  x1 <- sort(runif(n))
  x2 <- sort(runif(n))
  y <- 4 * (x1 - x2) + rnorm(n, sd = .5)

  # regress y on x1 and x2 under the shape-restriction: "decreasing-increasing"
  # with a penalty term = .1
  ans <- cgam(y ~ s.decr.incr(x1, x2), pen = .1)

# plot the constrained surface
  plotpersp(ans)
 }
}


\keyword{3D plot routine}
