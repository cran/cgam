\name{plotpersp}

\alias{plotpersp}

\title{Create a 3D Plot for a CGAM Object}

\description{Given an object of the cgam class, which has at least two non-parametrically modelled predictors, this routine will make a 3D plot of the fit with a set of two non-parametrically modelled predictors in the formula being the x and y labs. If there are more than two non-parametrically modelled predictors, any other such predictor will be evaluated at the largest value which is smaller than or equal to its median value.

If there is any categorical covariate and if the user specifies the argument categ to be a character representing a categorical covariate in the formula, then a 3D plot with multiple parallel surfaces, which represent the levels of a categorical covariate in an ascending order, will be created; otherwise, a 3D plot with only one surface will be created. Each level of a categorical covariate will be evaluated at its mode.

This routine is extended to make a 3D plot for an object fitted by warped-plane splines or triangle splines. Note that two non-parametrically modelled predictors specified in this routine must both be modelled as addtive components, or a pair of predictors forming an isotonic or convex surface without additivity assumption.

This routine is an extension of the generic R graphics routine persp. See the documentation below for more details. 

}

\usage{
  plotpersp(object,...)
}

\arguments{

\item{object}{An object of the cgam class with at least two non-parametrically modelled predictors.}
\item{...}{
  Additional arguments passed to the plotting method, such as \code{x1}, \code{x2}, and \code{control}. 
  These include:
  \describe{
    \item{\code{x1}}{Optional; either a character string or a vector matching the first nonparametric predictor. If omitted, the first predictor's name is used.}
    \item{\code{x2}}{Optional; either a character string or a vector matching the second nonparametric predictor. If omitted, the second predictor's name is used.}
    \item{\code{control}}{A list of graphical control parameters, typically created with \code{\link[cgam]{plotpersp_control}}. Controls aspects such as colors, axis limits, rotation angles, and grid resolution.}
  }
}
 
   %\item{...}{Additional arguments passed to \code{\link[graphics]{persp}}, such as `box`, `axes`, or `expand`.}
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

  ctl <- list(ngrid = 10, main = "Cgam Increasing Fit", 
  sub = "Categorical Variable: Sex", categ = "factor(sex)")
  plotpersp(fit, age, height, control = ctl)
  
  ctl <- list(ngrid = 10, main = "Cgam Smooth Increasing Fit", 
  sub = "Categorical Variable: Smoke", categ = "factor(smoke)")
  plotpersp(fit.s, age, height, control = ctl)

# Example 2.
  data(plasma)
  
  fit <- cgam(logplasma ~  s.decr(bmi) + s.decr(logdietfat) + s.decr(cholest) + s.incr(fiber) 
+ s.incr(betacaro) + s.incr(retinol) + factor(smoke) + factor(vituse), data = plasma) 

  ctl <- list(ngrid = 15, th = 120, ylab = "log(dietfat)",
            zlab = "est mean of log(plasma)", main = "Cgam Fit with the Plasma Data Set",
            sub = "Categorical Variable: Vitamin Use", categ = "factor(vituse)")
 
 plotpersp(fit, "bmi", "logdietfat", control = ctl)

# Example 3.
  data(plasma)
  addl <- 1:314*0 + 1 
  addl[runif(314) < .3] <- 2
  addl[runif(314) > .8] <- 4
  addl[runif(314) > .8] <- 3

  plasma$addl <- addl 
  ans <- cgam(logplasma ~ s.incr(betacaro, 5) + s.decr(bmi) +
                s.decr(logdietfat) + as.factor(addl), data = plasma)
  
  ctl <- list(th = 240, random = TRUE, categ = "as.factor(addl)")
  plotpersp(ans, "betacaro", "logdietfat", control = ctl)

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
