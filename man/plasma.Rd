\name{plasma}

\alias{plasma}

\docType{data}

\title{A Data Set for Cgam}

\description{This data set is used for the routine plotpersp. It contains 314 observations of blood plasma, beta carotene measurements along with several covariates. High levels of blood plasma and beta carotene are believed to be protective against cancer, and it is of interest to determine the relationships with covariates. }

\usage{data(plasma)}

\format{
  \describe{
    \item{\code{logplasma}}{A numeric vector of the logarithm of plasma levels.}
    \item{\code{betacaro}}{A numeric vector of dietary beta carotene consumed mcg per day.}
    \item{\code{bmi}}{A numeric vector of BMI values.}
    \item{\code{cholest}}{A numeric vector of cholesterol consumed mg per day.}
    \item{\code{dietfat}}{A numeric vector of the logarithm of grams of diet fat consumed per day.}
    \item{\code{fiber}}{A numeric vector of grams of fiber consumed per day.}
    \item{\code{retinol}}{A numeric vector of retinol consumed per day.}
    \item{\code{smoke}}{A numeric vector of smoking status (1=Never, 2=Former, 3=Current Smoker).}
    \item{\code{vituse}}{A numeric vector of vitamin use (1=Yes, fairly often, 2=Yes, not often, 3=No).}
  }
}

\source{\url{http://axon.cs.byu.edu/data/statlib/numeric/plasma_retinol.arff}}

\examples{data(plasma)}

\keyword{datasets}
