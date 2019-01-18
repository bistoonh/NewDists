\name{EEN}
\alias{dEEN}
\alias{pEEN}
\alias{qEEN}
\alias{rEEN}
\alias{mEEN}
\alias{vEEN}
\alias{sEEN}
\alias{kEEN}
\alias{entEEN}
\title{Extended Exp-Normal Distribution}
\description{
These functions provide information about the Extended Exp-Normal distribution.
dEEN gives the density,
pEEN gives the distribution function,
qEEN gives the quantile function,
rEEN generates random deviates,
mEEN gives the mean function,
vEEN gives the varince function,
sEEN gives the skewness function,
kEEN gives the kurtosis function, and
entEEN gives the Rrnyi or Shannon entropy function due to Alizadeh et al. (2018) specified by the pdf
\deqn{
f(x;\alpha ,\beta ,\mu ,\sigma )=\frac{\phi \left( \frac{x-\mu}{\sigma} \right)\Phi
{{\left( \frac{x-\mu}{\sigma} \right)}^{\alpha -1}}\left[ \alpha +\left( \beta -\alpha  \right)
\Phi {{\left( \frac{x-\mu}{\sigma}  \right)}^{\beta }} \right]}{{{\left[ \Phi {{\left( \frac{x-\mu}{\sigma}
 \right)}^{\alpha }}+1-\Phi {{\left(\frac{x-\mu}{\sigma}  \right)}^{\beta }} \right]}^{2}}},
}
where \eqn{\alpha, \beta, \sigma > 0}.}
\usage{
dEEN(x, alpha, beta, mu = 0, sigma = 1, log = FALSE)
pEEN(x, alpha, beta, mu = 0, sigma = 1)
qEEN(p, alpha, beta, mu = 0, sigma = 1)
rEEN(n, alpha, beta, mu = 0, sigma = 1)
mEEN(alpha, beta, mu = 0, sigma = 1)
vEEN(alpha, beta, mu = 0, sigma = 1)
sEEN(alpha, beta, mu = 0, sigma = 1)
kEEN(alpha, beta, mu = 0, sigma = 1)
entEEN(gamma, alpha, beta, mu = 0, sigma = 1, explain = FALSE)
}
\arguments{
  \item{x}{Scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{Scaler or vector of probabilities at which the quantile needs to be computed}
  \item{n}{Number of random numbers to be generated}
  \item{alpha}{The value of the first shape parameter. Must be positive and finite.}
  \item{beta}{The value of the second shape parameter. Must be positive and finite.}
  \item{mu}{Value of mean. Must be finite.}
  \item{sigma}{Value of standard deviations. Must be positive and finite.}
  \item{log}{Logical; if TRUE, probabilities p are given as log(p).}
  \item{gamma}{The gamma in Rrnyi entropy. if gamma = 1, the Shannon entropy is returned.}
  \item{explain}{Logical; if TRUE, explain Rrnyi or Shannon entropy is returned.}
  }

\value{
An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or
an object of the same length as \code{p}, giving the quantile values computed at \code{p} or
an object of the same length as \code{n}, giving the random numbers generated or
an object of the same length as \code{gamma}, giving the entropy (Rrnyi or Shannon) or
an object giving the values of mean, varince, skewness, or kurtosis.
 }
\references{
Alizadeh, Morad, Mahmoud Afshari, Bistoon Hosseini, and Thiago G. Ramires.
"Extended exp-G family of distributions: Properties, applications and simulation."
Communications in Statistics-Simulation and Computation (2018): 1-16.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
x = runif(10, min= -1 , max = 1)
dEEN(x, alpha = 0.5, beta = 1.2, mu = -2, sigma = 0.2, log = FALSE)
pEEN(x, alpha = 0.5, beta = 1.2, mu = 1, sigma = 1.1)
qEEN(x, alpha = 0.5, beta = 1.2, mu = .5, sigma = 11)
rEEN(n = 10, alpha = 0.5, beta = 1.2, mu = 0, sigma = 1)
mEEN(alpha = 0.5, beta = 1.2, mu = 10, sigma = .1)
vEEN(alpha = 0.5, beta = 1.2, mu = 2, sigma = 3)
sEEN(alpha = 0.5, beta = 1.2, mu = 1, sigma = 1)
kEEN(alpha = 0.5, beta = 1.2, mu = 0, sigma = 2)
entEEN(gamma = c(.5, 1, 1.5), alpha = 0.5, beta = 1.2, mu = 0, sigma = 1, explain = TRUE)

