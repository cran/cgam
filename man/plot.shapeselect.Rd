\name{plot.shapeselect}
\alias{plot.shapeselect}

\title{Plot Routine for a SHAPESELECT Fit}

\description{This routine is registered as the S3 plot method for a ShapeSelect fit.}

\usage{\method{plot}{shapeselect}(x, ...)}

\arguments{
  \item{x}{A fitted object of the ShapeSelect routine.}
  \item{...}{Furhter possible arguments.}
}


\value{The routine makes a plot for a ShapeSelect object, in which the xlab represents generations, i.e., iterations, a ShapeSelect routine goes through until it finds the best model by the genetic algorithm, and the ylab represents fitness values. Both the best fitness and the mean fitness in each generation is shown in the plot.}

\details{
Note that this routine only works on a ShapeSelect object fitted with the genetic algorithm; otherwise, if the object is found via brutal search, the routine will not work on it.
}

\author{Xiyue Liao}

\seealso{\code{\link{ShapeSelect}}}

\examples{
\dontrun{
  library(MASS)
  data(Rubber)

  # call ShapeSelect to find the best model by the genetic algorithm
  fit <- ShapeSelect(loss ~ shapes(hard, set = "s.9") + shapes(tens, set = "s.9"),
  data = Rubber, genetic = TRUE)

  # plot fitness values against generations
  plot(fit)

  # call plotpersp on the best fit 
  plotpersp(fit$best.fit, th = 120)
  
}
}

\keyword{plot routine}
