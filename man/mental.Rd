\name{mental}
\alias{mental}
\docType{data}
\title{Alachua County Study of Mental Impairment}

\description{The date set is from a study of mental health for a random sample of 40 adult residents of Alachua County, Florida. Mental impairment is an ordinal response with 4 categories: well, mild symptom formation, moderate symptom formation, and impaired, which are recorded as 1, 2, 3, and 4. Life event index is a composite measure of the number and severity of important life events that occurred with the past three years, e.g., birth of a child, new job, divorce, or death of a family member. It is an integer from 0 to 9. Another covariate is socio-economic status and it is measured as binary: high = 1, low = 0.  
}

\usage{data(mental)}

\format{
  \describe{
    \item{\code{mental}}{Mental impairment. It is an ordinal response with \eqn{4} categories recorded as \eqn{1}, \eqn{2}, \eqn{3}, and  \eqn{4}.}
    \item{\code{ses}}{Socio-economic status measured as binary: high = \eqn{1}, low = \eqn{0}.  
}
    \item{\code{life}}{Life event index. It is an integer from \eqn{0} to \eqn{9}.}
  }
}

\references{
		Agresti, A. (2010) Analysis of Ordinal Categorical Data, \emph{2nd ed. Hoboken, NJ: Wiley}.
}

\seealso{
  \code{\link{Ord}}
}

\examples{
	# proportional odds model example 
	data(mental)
	
	# model the relationship between the latent variable and life event index as increasing
	# socio-economic status is included as a binary covariate 
	fit.incr <- cgam(mental ~ incr(life) + ses, data = mental, family = Ord)

	# check the estimated probabilities P(mental = k), k = 1, 2, 3, 4
	probs.incr <- fitted(fit.incr)
	head(probs.incr)
}
\keyword{datasets}
