\name{GOGaU}
\alias{dGOGaU}
\alias{pGOGaU}
\alias{qGOGaU}
\alias{rGOGaU}
\alias{mGOGaU}
\alias{vGOGaU}
\alias{sGOGaU}
\alias{kGOGaU}
\alias{entGOGaU}
\title{Generalized Odd Gamma Uniform Distribution}
\description{
These functions provide information about the Generalized Odd Gamma Uniform distribution.
dGOGaU gives the density,
pGOGaU gives the distribution function,
qGOGaU gives the quantile function,
rGOGaU generates random deviates,
mGOGaU gives the mean function,
vGOGaU gives the variance function,
sGOGaU gives the skewness function,
kGOGaU gives the kurtosis function, and
entGOGaU gives the Rrnyi or Shannon entropy function due to Hosseini et al. (2018) specified by the pdf
\deqn{
f(x;\alpha ,\beta ,a,b) =\frac{\beta (b-a)^{\beta}(x-a)^{\alpha\beta-1}}{\Gamma(\alpha)
\left[ (b-a)^{\beta}-(x-a)^{\beta} \right]^{\alpha+1}} e^{\frac{-(x-a)^{\beta}}{(b-a)^{\beta}-(x-a)^{\beta}}}  ,a \le x \le b\,.
}
where \eqn{\alpha, \beta > 0} and  \eqn{a<b}.}
\usage{
dGOGaU(x, alpha, beta, a = 0, b = 1, log = FALSE)
pGOGaU(x, alpha, beta, a = 0, b = 1)
qGOGaU(p, alpha, beta, a = 0, b = 1)
rGOGaU(n, alpha, beta, a = 0, b = 1)
mGOGaU(alpha, beta, a = 0, b = 1)
vGOGaU(alpha, beta, a = 0, b = 1)
sGOGaU(alpha, beta, a = 0, b = 1)
kGOGaU(alpha, beta, a = 0, b = 1)
entGOGaU(gamma, alpha, beta, a = 0, b = 1, explain = FALSE)
}
\arguments{
  \item{x}{Scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{Scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{Number of random numbers to be generated}
  \item{alpha}{The value of the first shape parameter. Must be finite.}
  \item{beta}{The value of the second shape parameter. Must be finite.}
  \item{a, b}{Lower and upper limits of the distribution. Must be finite.}
  \item{log}{Logical; if TRUE, probabilities p are given as log(p).}
  \item{gamma}{The gamma in Rrnyi entropy. if gamma = 1, the Shannon entropy is returned.}
  \item{explain}{Logical; if TRUE, explain Rrnyi or Shannon entropy is returned.}
  }

\value{
An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or
an object of the same length as \code{p}, giving the quantile values computed at \code{p} or
an object of the same length as \code{n}, giving the random numbers generated or
an object of the same length as \code{gamma}, giving the entropy (Rrnyi or Shannon) or
an object giving the values of mean, variance, skewness, or kurtosis.
 }
\references{Hosseini, Bistoon, Mahmoud Afshari, and Morad Alizadeh. "The Generalized Odd Gamma-G Family of Distributions: Properties and Applications." Austrian Journal of Statistics 47.2 (2018): 69-89.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{x=runif(10,min=0,max=1)
dGOGaU(x, alpha = 0.5, beta = 1.2, a = 0, b = 1, log = FALSE)
pGOGaU(x, alpha = 0.5, beta = 1.2, a = 0, b = 1)
qGOGaU(x, alpha = 0.5, beta = 1.2, a = 0, b = 1)
rGOGaU(n = 10, alpha = 0.5, beta = 1.2, a = 0, b = 1)
mGOGaU(alpha = 0.5, beta = 1.2, a = 0, b = 1)
vGOGaU(alpha = 0.5, beta = 1.2, a = 0, b = 1)
sGOGaU(alpha = 0.5, beta = 1.2, a = 0, b = 1)
kGOGaU(alpha = 0.5, beta = 1.2, a = 0, b = 1)
entGOGaU(gamma = c(.5, 1, 1.5), alpha = 0.5, beta = 1.2, a = 0, b = 1, explain = TRUE)

