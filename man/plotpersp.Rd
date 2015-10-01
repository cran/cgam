\name{plotpersp}

\alias{plotpersp}

\title{Create a 3D Plot for a Cgam Object}

\description{Given an object of the cgam class, which has at least two non-parametrically modelled predictors, this routine will make a 3D plot of the fit with a set of two non-parametrically modelled predictors in the formula being the x and y labs. If there are more than two non-parametrically modelled predictors, any other such predictor will be evaluated at the largest value which is smaller than or equal to its median value.

If there is any categorical covariate and if the user specifies the argument categ to be a character representing a categorical covariate in the formula, then a 3D plot with multiple parallel surfaces, which represent the levels of a categorical covariate in an ascending order, will be created; otherwise, a 3D plot with only one surface will be created. Each level of a categorical covariate will be evaluated at its mode.

This routine is an extension of the generic R graphics routine persp.}

\usage{plotpersp(object, x1, x2, data = NULL, surface = "mu", categ = NULL, 
col = NULL, random = FALSE, x_grid = 20, y_grid = 20, at = "median", 
xlim = range(x1), ylim = range(x2), zlim = NULL, xlab = NULL, ylab = NULL, 
zlab = NULL, main = NULL, sub = NULL, th = NULL, phi = 15, r = sqrt(3), 
d = 1, scale = TRUE, expand = 1, border = NULL, ltheta = NULL, lphi = 0, 
shade = NA, box = TRUE, axes = TRUE, nticks = 5, ticktype = "detailed")}

\arguments{
  \item{object}{An object of the cgam class with at least two non-parametrically modelled predictors.}
  \item{x1}{A non-parametrically modelled predictor in a cgam fit.}
  \item{x2}{A non-parametrically modelled predictor in a cgam fit.}
  \item{data}{The data frame based on which the user get a cgam fit.}
  \item{surface}{The type of the surface of a 3D plot. If surface == "mu", then the surface of the estimated mean value of a cgam fit will be plotted; if surface == "eta", then the surface of the estimated systematic component value of a cgam fit will be plotted. The default is surface = "mu".}
 \item{categ}{Optional categorical covariate(s) in a cgam fit. If there is any categorical covariate and if the user specifies the argument categ to be a character representing a categorical covariate in the formula, then a 3D plot with multiple parallel surfaces, which represent the levels of a categorical covariate in an ascending order, will be created; otherwise, a 3D plot with only one surface will be created. Each level of a categorical covariate will be evaluated at its mode. The default is categ = NULL.}
  \item{col}{The color(s) of a 3D plot created by plotpersp. If col == NULL, "white" will be used when there is only one surface in the plot, and a sequence of colors will be used in a fixed order when there are multiple parallel surfaces in the plot. For example, when there are two surfaces, the lower surface will be in the color "peachpuff", and the higher surface will be in the color "lightblue". The default is col =  NULL.}
  \item{random}{A logical scalar. If random == TRUE, color(s) for a 3D plot will be randomly chosen from ten colors, namely, "peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke"; otherwise, "white" will be used when there is only one surface in the plot, and a sequence of colors will be used in a fixed order when there are multiple parallel surfaces in the plot.}
  \item{x_grid}{This is a positive integer specifying how dense the \eqn{x} grid will be. The default is x_grid = \eqn{20}.}
  \item{y_grid}{This is a positive integer specifying how dense the \eqn{y} grid will be. The default is y_grid = \eqn{20}.}
  \item{at}{When a user specifies a set of two non-parametrically modelled predictors to make a 3D plot, any other non-parametrically modelled predictor in the cgam formula will be evaluated at the largest value which is smaller than or equal to its median value.}
  \item{xlim}{The xlim argument inherited from the persp routine.}
  \item{ylim}{The ylim argument inherited from the persp routine.}
  \item{zlim}{The zlim argument inherited from the persp routine.}
  \item{xlab}{The xlab argument inherited from the persp routine.}
  \item{ylab}{The ylab argument inherited from the persp routine.}
  \item{zlab}{The zlab argument inherited from the persp routine.}
  \item{main}{The main argument inherited from the persp routine.}
  \item{sub}{The sub argument inherited from the persp routine.}
  \item{th}{The theta argument inherited from the persp routine.}
  \item{phi}{The phi argument inherited from the persp routine.}
  \item{r}{The r argument inherited from the persp routine.}
  \item{d}{The d argument inherited from the persp routine.}
  \item{scale}{The scale argument inherited from the persp routine.}
  \item{expand}{The expand argument inherited from the persp routine.}
  \item{border}{The border argument inherited from the persp routine.}
  \item{ltheta}{The ltheta argument inherited from the persp routine.}
  \item{lphi}{The lphi argument inherited from the persp routine.}
  \item{shade}{The shade argument inherited from the persp routine.}
  \item{box}{The box argument inherited from the persp routine.}
  \item{axes}{The axes argument inherited from the persp routine.}
  \item{nticks}{The nticks argument inherited from the persp routine.}
  \item{ticktype}{The ticktype argument inherited from the persp routine.}
}

\value{The routine plotpersp returns a 3D plot of an object of the cgam class. The \eqn{x} lab and \eqn{y} lab represent a set of non-parametrically modelled predictors used in a cgam formula. The \eqn{z} lab represents the estimated mean value or the estimated systematic component value.}

\references{The official documentation for the generic R routine persp: \url{http://stat.ethz.ch/R-manual/R-patched/library/graphics/html/persp.html}}

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

  fit = cgam(y ~ incr(age) + incr(height) + sex + smoke, nsim = 0)
  fit1 = cgam(y ~ s.incr(age) + s.incr(height) + sex + smoke, nsim = 0)

  plotpersp(fit, age, height, x_grid = 10, y_grid = 10, main = "Cgam Increasing Fit", 
  sub = "Categorical Variable: Sex", categ = "sex")
  plotpersp(fit1, age, height, x_grid = 10, y_grid = 10, main = "Cgam Smooth Increasing Fit", 
  sub = "Categorical Variable: Smoke", categ = "smoke")

# Example 2.
  data(plasma)

  # extract the variables
  y <- plasma$logplasma
  bmi <- plasma$bmi
  dietfat <- plasma$dietfat
  cholest <- plasma$cholest
  fiber <- plasma$fiber
  betacaro <- plasma$betacaro
  retinol <- plasma$retinol
  smoke <- plasma$smoke
  vituse <- plasma$vituse

  fit <- cgam(y ~  s.decr(bmi) + s.decr(dietfat) + s.decr(cholest) + s.incr(fiber) 
+ s.incr(betacaro) + s.incr(retinol) + smoke + vituse) 

  plotpersp(fit, bmi, dietfat, x_grid = 10, y_grid = 10, th = 120, ylab = "log(dietfat)", 
zlab = "est mean of log(plasma)", main = "Cgam Fit with the Plasma Data Set", 
sub = "Categorical Variable: Vitamin Use", categ = "vituse")

# Example 3.
  data(plasma)
  addl <- 1:314*0 + 1 
  addl[runif(314) < .3] <- 2
  addl[runif(314) > .8] <- 4
  addl[runif(314) > .8] <- 3

  ans <- cgam(logplasma ~ s.incr(betacaro, 5) + s.decr(bmi) + s.decr(dietfat) 
+ as.factor(addl), data = plasma)
  plotpersp(ans, betacaro, dietfat, th = 240, random = TRUE, 
categ = "as.factor(addl)", data = plasma)
}

\keyword{3D plot routine}
