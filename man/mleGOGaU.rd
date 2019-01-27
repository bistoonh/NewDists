\name{mleGOGaU}
\alias{mleGOGaU}

\title{The MLEs of GOGaU}
\description{
This function calculates the MLEs of GOGaU distribution.
}

\usage{
mleGOGaU(x, par0 = c(1 , 1, min(x)-1, max(x)+1), a = NA, b = NA, fitplot = TRUE)
}
\arguments{
  \item{x}{The value of the first shape parameter. Must be finite.}
  \item{par0}{Initial values for the parameters to be optimized over.}
  \item{a, b}{Values for fixed lower or upper limits of the GOGaU distribution.
  If \code{a}/\code{b} is equal to \pkg{NA}, it is estimated as the unknown parameter.}
  
  \item{fitplot}{logical; if TRUE, histogram and fitted GOGaU density is drawn.}
	}
\value{
  \item{par}{The MLEs of parameters.}
  \item{loglike}{The log-likelihood value corresponding to MLEs.}
  \item{convergence}{An integer code. 0 indicates successful completion.
					To see possible error codes, see the \link{optim} function.
					}
  \item{W}{The statistic Cramér-von Misses.}
  \item{A}{The statistic Anderson Darling.}
  \item{KS}{Kolmogorov Smirnov test.}
  \item{AIC}{The value of Akaike Information Criterion.}
  \item{BIC}{The value of Bayesian Information Criterion.}
 }
\references{
Bistoon Hosseini, Mahmoud Afshari, and Morad Alizadeh. "The Generalized Odd Gamma-G Family of Distributions:
Properties and Applications." Austrian Journal of Statistics 47.2 (2018): 69-89.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
x = rGOGaU(n = 1000, alpha = 0.5, beta = 1.5, a = 0, b = 4.5)
mleGOGaU(x = x, a = NA, b = 4.5)
}

