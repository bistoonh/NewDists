\name{mleEEN}
\alias{mleEEN}

\title{The MLEs of EEN}
\description{
This function calculates the MLEs of EEN distribution.
}

\usage{
mleEEN(x, par0 = c(1, 1, mean(x), sd(x)), fitplot = TRUE)
}
\arguments{
  \item{x}{The value of the first shape parameter. Must be finite.}
  \item{par0}{Initial values for the parameters to be optimized over.}
  If \code{a}/\code{b} is equal to \pkg{NA}, it is estimated as the unknown parameter.}
  }
  \item{fitplot}{Logical; if TRUE, histogram and fitted EEN density is drawn.}

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
Alizadeh, Morad, Mahmoud Afshari, Bistoon Hosseini, and Thiago G. Ramires.
"Extended exp-G family of distributions: Properties, applications and simulation."
Communications in Statistics-Simulation and Computation (2018): 1-16.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
x = rEEN(n = 1000, alpha = 1.5, beta = 1.5, mu = 0, sigma = 1)
mleEEN(x)
}

