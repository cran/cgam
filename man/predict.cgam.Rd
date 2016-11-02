\name{predict.cgam}
\alias{predict.cgam}
\title{Prediction Routine for a CGAM Fit}

\description{This routine is registered as the S3 prediction method for a cgam fit.}

\usage{\method{predict}{cgam}(object, newData,...)}

\arguments{
  \item{object}{A fitted object of the cgam routine.}
  \item{newData}{A data frame containing the values of the predictor(s), at which the prediction is made, in the formula of a cgam fit.}
  \item{...}{Furhter possible arguments.}
}

\details{
  The prediction routine mimics predict.lm and predict.glm. To use this routine, the user should define ``newData'' as a data frame, in which the name of each column should match a name appearing in the formula of a cgam fit.

  Note that for now, we only have a prediction routine for a cgam fit, not a wps fit. 
}

\value{
  \item{muhat}{A vector of predicted values for the predictor(s) defined in the newData argument.}
}

\references{
  Meyer, M. C. (2013a) Semi-parametric additive constrained regression. \emph{Journal of Nonparametric Statistics \bold{25(3)}}, 715

  Meyer, M. C. (2008) Inference using shape-restricted regression splines. \emph{Annals of Applied Statistics \bold{2(3)}}, 1013--1033.
}

\author{Xiyue Liao and Mary C. Meyer}

\seealso{\code{\link{cgam}}}

\examples{
# Example 1.
  library(SemiPar)
  data(bpd)

  # an indicator of presence of bronchopulmonary dysplasia (BPD): 0 = absent, 1 = present
  BPD <- bpd$BPD
 
  # birthweight of baby (grammes)
  birthweight <- bpd$birthweight
  fit <- cgam(BPD ~ s.decr(birthweight, space = "Q"), family = binomial(), nsim = 0)

  # make a data frame which is the birthweight mean and call the prediction routine
  xp = mean(birthweight)
  newd = data.frame(birthweight = xp)
  mp = predict(fit, newd)

  # check the prediction in the plot
  plot(birthweight, BPD, type = "n", xlab = "Birthweight (grams)", 
  ylab = "Bronchopulmonary dysplasia", cex.lab = 1.4)
  rug(birthweight[BPD == 0])
  rug(birthweight[BPD == 1], side = 3)
  ord = order(birthweight)
  lines(sort(birthweight), fit$muhat[ord], lwd = 2)
  points(xp, mp, col = 2, pch = 16, cex = 1.3)
  abline(v = xp, lty = 2)
  legend("right", bty = "n", c("predicted value for the birthweight mean"), 
  col = 2, pch = 16, cex = .6)

# Example 2.
  # an example with a categorical covariate
  data(feet)
  l <- feet$length 
  w <- feet$width
  sex <- feet$sex
  fit <- cgam(w ~ s.incr(l) + factor(sex), nsim = 0)

  # use the first five 'l' values and five new 'sex' values to make prediction based on 'fit'
  newd <- data.frame(l = l[1:5], sex = c("B", "G", "G", "B", "G"))
  mp <- predict(fit, newd)

  # check the prediction
  mp

  # check the first five 'sex' values in the 'feet' data set
  sex[1:5]

  # if we adjust the prediction with the estimated coefficient for sex = 'G', 
  # then the prediction should be the same as the fit for the first five observations
  z_add <- fit$zcoefs[2] 
  mp[c(2, 3, 5)] <- mp[c(2, 3, 5)] - z_add

  # check if they are equal
  fit$muhat[1:5]
  #[1] 9.001077 9.260961 9.025589 9.208811 9.182691
  mp
  #[1] 9.001077 9.260961 9.025589 9.208811 9.182691
}

\keyword{prediction routine}
